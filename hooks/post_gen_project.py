#!/usr/bin/env python
import os
import shutil
import sys

from pm4ngs.utils import copy_rawdata_to_project

if __name__ == '__main__':
    DATASET = '{{ cookiecutter.dataset_name }}'
    DATASET_DIR = os.path.join(os.path.realpath(os.path.curdir), 'data', DATASET)
    SAMPLE_TABLE_FILE = os.environ.get('PM4NGS_SAMPLE_TABLE', None)
    COPY_RAWDATA = os.environ.get('PM4NGS_COPY_RAWDATA', None)

    if SAMPLE_TABLE_FILE and COPY_RAWDATA:
        print('Copying file {}  to {}'.format(
            SAMPLE_TABLE_FILE, os.path.join(DATASET_DIR, 'sample_table.csv')
        ))
        shutil.copyfile(SAMPLE_TABLE_FILE, os.path.join(DATASET_DIR, 'sample_table.csv'))
        copy_rawdata_to_project(COPY_RAWDATA, DATASET_DIR)
        print(' Done')
    else:
        print('Error reading user env')
        print('PM4NGS_SAMPLE_TABLE: ' + str(SAMPLE_TABLE_FILE))
        print('PM4NGS_COPY_RAWDATA: ' + str(COPY_RAWDATA))
        sys.exit(-1)
