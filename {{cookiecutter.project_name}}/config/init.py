import os
import re
import json
import pandas
import time
import math
import pickle
import zipfile
import uuid
import gzip
import distutils.spawn
import numpy as np
import scipy.stats as stats
import seaborn as sns
import locale

from multiprocessing import Pool

from Bio import SeqIO
from datetime import datetime
import networkx as nx

from pm4ngs.jupyterngsplugin.files.fastq.fastqc import parse_fastqc_zip
from pm4ngs.jupyterngsplugin.markdown.utils import find_file_print_link_size
from pm4ngs.jupyterngsplugin.utils.errors import check_cwl_command_log
from pm4ngs.jupyterngsplugin.utils.run_command import run_command
from pm4ngs.jupyterngsplugin.utils.working_dir import working_dir
from pm4ngs.jupyterngsplugin.utils.yaml_utils import write_to_yaml
from pm4ngs.jupyterngsplugin.utils.yaml_utils import load_from_yaml

from pm4ngs.jupyterngsplugin.markdown.utils import find_file_print_link_size, get_link_size
from pm4ngs.jupyterngsplugin.utils.load_content_dict import load_content_dict_line

from IPython.display import display

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as lines
import matplotlib.cm as cm
from matplotlib.gridspec import GridSpec

from matplotlib import font_manager as fm, rcParams

from IPython.display import HTML
from IPython.display import display, Markdown, Latex

###############################################################
#
#    Project global paths
#
###############################################################

WORKDIR = os.path.dirname(os.path.abspath(os.path.dirname(__file__)))
CONFIG = os.path.join(WORKDIR,'config')
DATA = os.path.join(WORKDIR,'data')
BIN = os.path.join(WORKDIR,'bin')
RESULTS = os.path.join(WORKDIR,'results')
NOTEBOOKS = os.path.join(WORKDIR,'notebooks')
SRC = os.path.join(WORKDIR,'src')
TMP = os.path.join(WORKDIR,'tmp')


###############################################################
#
#    Dataset (experiment) to analyze
#
# The path is $WORKDIR/data/$DATASET
#
# To use multiple datasets (experiments) this variable should be overwritten
# in the notebooks
#
###############################################################

DATASET = '{{ cookiecutter.dataset_name }}'

###############################################################
#
#    Checking installed tools
#
###############################################################

if not distutils.spawn.find_executable('gcloud'):
    raise ImportError('gcloud not in path.\nInstall: google-cloud-sdk\n')

if not distutils.spawn.find_executable('elastic-blast.py'):
    raise ImportError('elastic-blast not in path.\nInstall: elastic-blast\n')

if not distutils.spawn.find_executable('kubectl'):
    raise ImportError('kubectl not in path.\nInstall: kubernetes-client, version 1.18.8\n')

###############################################################
#
#    Utils functions
#
###############################################################


def getActionTime(data, actionId):
    ts = None
    te = None
    for l in data:
        if 'containerStarted' in l:
            if l['containerStarted']['actionId'] == actionId:
                ts = datetime.strptime(l['timestamp'].split('.')[0], "%Y-%m-%dT%H:%M:%S")
        if 'containerStopped' in l:
            if l['containerStopped']['actionId'] == actionId:
                te = datetime.strptime(l['timestamp'].split('.')[0], "%Y-%m-%dT%H:%M:%S")
        if ts and te:
            return te - ts
    return None


def get_gpc_starttimestamp(logs):
    ts = ""
    for e in logs['metadata']['events']:
        if 'workerAssigned' in e and 'machineType' in e['workerAssigned']:
            ts = e['timestamp']
            break
    return ts


def parse_gcp_json(logs, f, blastdb_action, cwl_action):
    if 'done' in logs and logs['done'] == True:
        ts = get_gpc_starttimestamp(logs)
        ts = datetime.strptime(ts.split('.')[0], "%Y-%m-%dT%H:%M:%S")
        te = datetime.strptime(logs['metadata']['endTime'].split('.')[0], "%Y-%m-%dT%H:%M:%S")
        elapsed = te - ts
        blastdb = getActionTime(logs['metadata']['events'], blastdb_action)
        cwl = getActionTime(logs['metadata']['events'], cwl_action)
        return [f, elapsed, blastdb, cwl]
    return []


def successors(g, O):
    """
    Extract ancestors nodes from an starting node
    :param g: starting node name
    :param O: Graph
    :return: a set with node names
    """
    result = {g}
    for o in O.successors(g):
        result.update(successors(o, O))
    return result


def findNode(O, id):
    nodes = [y for x,y in O.nodes(data=True) if y['id']==id]
    if nodes:
        a = ""
        for i in nx.shortest_path(O, source="1", target=id)[2:]:  
            ns = [y for x,y in O.nodes(data=True) if y['id']==i]
            if ns:
                if a:
                    a += "; "
                a += ns[0]['name']
        return nodes[0], a
    return None, None


def parse_nodes_file(node_file, taxid):
    with open(node_file, 'r') as fin:
        for line in fin:
            f = line.strip().split('\t|\t')
            node = {}
            node['id'] = f[0]
            node['name'] = taxid[f[0]]['scientific name']
            edge = ()
            if f[1] != node['id']:
                edge = (f[1], node['id'])
            yield (f[0], node), edge


def parse_tax_name_file(name_file):
    tax_id = {}
    with open(name_file, 'r') as fin:
        for line in fin:
            line = line.strip()
            f = line.split('\t|\t')
            if f[0] not in tax_id:
                tax_id[f[0]] = {}
            tax_id[f[0]][f[3].replace('\t|','').strip()] = f[1]
    return tax_id
