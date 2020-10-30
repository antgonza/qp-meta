# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# This file contains functions related to generating default parameters for
# meta
# ------------------------------------------------------------------------------

import os
from os.path import join, isdir
import pandas as pd
from biom import Table

ALIGNERS = [
    "bowtie2"]


def get_dbs(db_folder):
    dbs = {}
    # Loop through the databases and create a dict of them
    for folder in os.listdir(db_folder):
        folder_path = join(db_folder, folder)
        if isdir(folder_path):
            dbs[folder] = folder_path

    return(dbs)


def get_dbs_list(db_folder):
    dbs = []
    # Loop through the databases and create a list string
    for folder in sorted(os.listdir(db_folder)):
        folder_path = join(db_folder, folder)
        if isdir(folder_path):
            dbs.append(folder_path)
    dbs_formatted = (', '.join('"' + item + '"' for item in dbs))

    return(dbs_formatted)


def generate_meta_dflt_params():
    dflt_param_set = {}
    db_parent_path = os.environ["QC_meta_DB_DP"]
    # Get a the databases available and the database name
    dbs = get_dbs(db_parent_path)
    # Create dict with command options per database
    for db in dbs:
        for aligner in ALIGNERS:
            dflt_param_set[db+'_'+aligner] = {'Database': dbs[db],
                                              'Aligner tool': aligner,
                                              'Percent identity': 0.95,
                                              'Capitalist': False,
                                              'Number of threads': 15}

    return(dflt_param_set)


def readfq(fp):  # this is a generator function
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for line in fp:  # search for the start of the next record
                if line[0] in '>@':  # fasta/q header line
                    last = line[:-1]  # save this line
                    break
        if not last:
            break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for line in fp:  # read the sequence
            if line[0] in '@+>':
                last = line[:-1]
                break
            seqs.append(line[:-1])
        if not last or last[0] != '+':  # this is a fasta record
            yield name, ''.join(seqs), None  # yield a fasta record
            if not last:
                break
        else:  # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for line in fp:  # read the quality
                seqs.append(line[:-1])
                leng += len(line) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs)  # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield name, seq, None  # yield a fasta record instead
                break


def meta_db_functional_parser(db_path):
    # Metadata file path
    md_fp = join(db_path, 'metadata.yaml')
    metadata = pd.read_csv(md_fp, sep=':', index_col=0)
    func_prefix = metadata.loc['function'].values[0].strip()
    fp_array = {
        'enzyme': join(db_path, '%s-enzyme-annotations.txt' % func_prefix),
        'module': join(db_path, '%s-module-annotations.txt' % func_prefix),
        'pathway': join(db_path, '%s-pathway-annotations.txt' % func_prefix)}

    return fp_array


def meta_parse_enzyme_table(f):
    md = pd.read_csv(
        f, sep='\t', header=None, error_bad_lines=False, warn_bad_lines=False)
    md.set_index(0, inplace=True)
    metadata = {}
    for i, row in md.iterrows():
        metadata[i] = {'taxonomy': [x for x in row.values]}
    return(metadata)


def meta_parse_module_table(f):
    md = pd.read_csv(
        f, sep='\t', header=None, error_bad_lines=False, warn_bad_lines=False)
    metadata = {}
    for i, row in md.iterrows():
        module = row[4].split('  ')[0]
        name = row[4].split('  ')[1]
        if module not in metadata:
            metadata[module] = {'taxonomy': [row[1], row[2], row[3], name]}
    return(metadata)


def meta_parse_pathway_table(f):
    md = pd.read_csv(
        f, sep='\t', header=None, error_bad_lines=False, warn_bad_lines=False)
    metadata = {}
    for i, row in md.iterrows():
        pathway = row[4]
        if pathway not in metadata:
            metadata[pathway] = {'taxonomy': [row[1], row[2], row[3]]}
    return(metadata)


def import_meta_biom(f, annotation_table=None,
                     annotation_type=None, names_to_taxonomy=False):
    import_funcs = {'module': meta_parse_module_table,
                    'pathway': meta_parse_pathway_table,
                    'enzyme': meta_parse_enzyme_table}

    table = pd.read_csv(f, sep='\t', index_col=0)

    bt = Table(table.values,
               observation_ids=list(map(str, table.index)),
               sample_ids=list(map(str, table.columns)))

    if names_to_taxonomy:
        metadata = {
            x: {'taxonomy': x.split(';')} for x in bt.ids(axis='observation')}
        bt.add_metadata(metadata, axis='observation')

    if annotation_table is not None:
        metadata = import_funcs[annotation_type](annotation_table)
        bt.add_metadata(metadata, axis='observation')

    return(bt)
