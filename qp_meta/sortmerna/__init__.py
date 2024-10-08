# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from qiita_client import QiitaCommand

from .sortmerna import sortmerna

__all__ = ['sortmerna']

# Defining the Sortmerna command
req_params = {'input': ('artifact', ['per_sample_FASTQ'])}
opt_params = {}
outputs = {'Non-ribosomal reads': 'per_sample_FASTQ',
           'Ribosomal reads': 'per_sample_FASTQ'}

# defining default parameter set AKA what's going to be shown to the user
# as options for the command
dflt_param_set = {
    'Ribosomal read filtering': {}
}

sortmerna_cmd = QiitaCommand(
    'SortMeRNA v4.3.7', "Ribosomal read filtering", sortmerna,
    req_params, opt_params, outputs, dflt_param_set)
