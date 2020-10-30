# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from qiita_client import QiitaCommand
from .meta import meta
from .utils import (generate_meta_dflt_params, get_dbs_list)
from os import environ


__all__ = ['meta']

# Define the meta command
default_db_list = get_dbs_list(environ["QC_meta_DB_DP"])
req_params = {'input': ('artifact', ['per_sample_FASTQ'])}
opt_params = {
    # database
    'Database': ["choice: [%s]" % default_db_list,
                 # making the first option default and rm quotes
                 default_db_list.split(',')[0].strip('"')],
    # aligner
    'Aligner tool': ['choice:[' +
                     '"bowtie2"]', 'bowtie2'],
    # threads
    'Number of threads': ['integer', '10'],
    'Capitalist': ['boolean', 'False'],
    'Percent identity': ['float', '0.95'],
    }
outputs = {
    'meta Alignment Profile': 'BIOM',
    'Taxonomic Predictions - phylum': 'BIOM',
    'Taxonomic Predictions - genus': 'BIOM',
    'Taxonomic Predictions - species': 'BIOM',
    'meta - per genome': 'BIOM',
    'meta - per gene': 'BIOM',
    }
dflt_param_set = generate_meta_dflt_params()

meta_cmd = QiitaCommand(
    'meta v0.1.1', "Functional and Taxonomic Predictions", meta,
    req_params, opt_params, outputs, dflt_param_set)
