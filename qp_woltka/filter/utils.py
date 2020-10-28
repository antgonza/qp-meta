# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# This file contains functions related to generating default parameters for
# QC_Filter
# ------------------------------------------------------------------------------

import os
from os.path import join, isdir


def get_dbs(db_folder):
    dbs = {}
    # Loop through the databases and create a dict of them
    for folder in os.listdir(db_folder):
        folder_path = join(db_folder, folder)
        if isdir(folder_path):
            dbs[folder] = join(folder_path, folder)

    return(dbs)


def get_dbs_list(db_folder):
    dbs = []
    # Loop through the databases and create a list string
    for folder in sorted(os.listdir(db_folder)):
        folder_path = join(db_folder, folder)
        if isdir(folder_path):
            dbs.append(join(folder_path, folder))
    dbs_formatted = (', '.join('"' + item + '"' for item in dbs))

    return(dbs_formatted)


def generate_filter_dflt_params():
    dflt_param_set = {}
    db_parent_path = os.environ["QC_FILTER_DB_DP"]
    # Get a the databases available and the database name
    dbs = get_dbs(db_parent_path)
    # Create dict with command options per database
    for db in dbs:
        dflt_param_set[db] = {'Bowtie2 database to filter': dbs[db],
                              'Number of threads': 15}

    return(dflt_param_set)
