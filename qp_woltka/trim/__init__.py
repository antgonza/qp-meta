# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from qiita_client import QiitaCommand

from .trim import trim

__all__ = ['trim']

# Define the trim command
req_params = {'input': ('artifact', ['per_sample_FASTQ'])}
opt_params = {
    # 3' adapter
    'Fwd read adapter': ['string', 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'],
    # 3' adapter for rev
    'Rev read adapter': ['string', 'GATCGGAAGAGCGTCGTGTAGGGAAAGGAGTGT'],
    # 3' quality cutoff
    'Trim low-quality bases': ['integer', '15'],
    # min length after trimming
    'Minimum trimmed read length': ['integer', '80'],
    # drop pairs whose mates are filtered out
    'Pair-end read required to match': ['choice:["any", "both"]', 'any'],
    # maximum Ns to drop sequence
    'Maximum number of N bases in a read to keep it': ['integer', '80'],
    # trim Ns on end of read
    'Trim Ns on ends of reads': ['boolean', True],
    # Threads used
    'Number of threads used': ['integer', '15'],
    # NextSeq-specific quality trimming
    'NextSeq-specific quality trimming': ['boolean', False],
    }
outputs = {'Adapter trimmed files': 'per_sample_FASTQ'}
dflt_param_set = {
    'KAPA HyperPlus with iTru': {
        'Fwd read adapter': 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC',
        'Rev read adapter': 'GATCGGAAGAGCGTCGTGTAGGGAAAGGAGTGT',
        'Trim low-quality bases': 15,
        'Minimum trimmed read length': 80,
        'Pair-end read required to match': 'any',
        'Maximum number of N bases in a read to keep it': 80,
        'Trim Ns on ends of reads': True,
        'NextSeq-specific quality trimming': False,
        'Number of threads used': 15
        }
}

trim_cmd = QiitaCommand(
    'Atropos v1.1.24', "Sequence QC - adapter trimming", trim,
    req_params, opt_params, outputs, dflt_param_set)
