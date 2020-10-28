# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from os.path import join
from qp_woltka.utils import (
    _format_params, make_read_pairs_per_sample,
    _run_commands, _per_sample_ainfo)

ATROPOS_PARAMS = {
    'adapter': 'Fwd read adapter', 'A': 'Rev read adapter',
    'quality-cutoff': 'Trim low-quality bases',
    'minimum-length': 'Minimum trimmed read length',
    'pair-filter': 'Pair-end read required to match',
    'max-n': 'Maximum number of N bases in a read to keep it',
    'trim-n': 'Trim Ns on ends of reads', 'threads': 'Number of threads used',
    'nextseq-trim': 'NextSeq-specific quality trimming'}


def generate_trim_commands(forward_seqs, reverse_seqs, map_file,
                           out_dir, parameters):
    """Generates the QC_Trim commands

    Parameters
    ----------
    forward_seqs : list of str
        The list of forward seqs filepaths
    reverse_seqs : list of str
        The list of reverse seqs filepaths
    map_file : str
        The path to the mapping file
    out_dir : str
        The job output directory
    parameters : dict
        The command's parameters, keyed by parameter name

    Returns
    -------
    cmds: list of str
        The QC_Trim commands
    samples: list of tup
        list of 4-tuples with run prefix, sample name, fwd read fp, rev read fp

    Notes
    -----
    Currently this is requiring matched pairs in the make_read_pairs_per_sample
    step but implicitly allowing empty reverse reads in the actual command
    generation. This behavior may allow support of situations with empty
    reverse reads in some samples, for example after trimming and QC.
    """
    # we match filenames, samples, and run prefixes
    samples = make_read_pairs_per_sample(forward_seqs, reverse_seqs, map_file)
    cmds = []

    param_string = _format_params(parameters, ATROPOS_PARAMS)

    for run_prefix, sample, f_fp, r_fp in samples:
        if r_fp is None:
            cmds.append("atropos trim %s -o %s  -se %s" % (
                param_string, join(out_dir, '%s.R1.fastq.gz' % run_prefix),
                f_fp))
        else:
            cmds.append('atropos trim %s -o %s -p %s -pe1 %s -pe2 %s'
                        % (param_string, join(out_dir, '%s.R1.fastq.gz' %
                           run_prefix), join(out_dir, '%s.R2.fastq.gz' %
                           run_prefix), f_fp, r_fp))
    return cmds, samples


def trim(qclient, job_id, parameters, out_dir):
    """Run Atropos with the given parameters

    Parameters
    ----------
    qclient : tgp.qiita_client.QiitaClient
        The Qiita server client
    job_id : str
        The job id
    parameters : dict
        The parameter values to run split libraries
    out_dir : str
        The path to the job's output directory

    Returns
    -------
    bool, list, str
        The results of the job
    """
    # Step 1 get the rest of the information need to run Atropos
    qclient.update_job_step(job_id, "Step 1 of 4: Collecting information")
    artifact_id = parameters['input']
    del parameters['input']

    # Get the artifact filepath information
    artifact_info = qclient.get("/qiita_db/artifacts/%s/" % artifact_id)
    fps = artifact_info['files']

    # Get the artifact metadata
    prep_info = qclient.get('/qiita_db/prep_template/%s/'
                            % artifact_info['prep_information'][0])
    qiime_map = prep_info['qiime-map']

    # Step 2 generating command atropos
    qclient.update_job_step(job_id, "Step 2 of 4: Generating"
                                    " QC_Trim commands")
    rs = fps['raw_reverse_seqs'] if 'raw_reverse_seqs' in fps else []
    commands, samples = generate_trim_commands(fps['raw_forward_seqs'],
                                               rs, qiime_map, out_dir,
                                               parameters)

    # Step 3 execute atropos
    len_cmd = len(commands)
    msg = "Step 3 of 4: Executing QC_Trim job (%d/{0})".format(len_cmd)
    success, msg = _run_commands(qclient, job_id, commands, msg, 'QC_Trim')
    if not success:
        return False, None, msg

    # Step 4 generating artifacts
    msg = "Step 4 of 4: Generating new artifacts (%d/{0})".format(len_cmd)
    suffixes = ['%s.R1.fastq.gz', '%s.R2.fastq.gz']
    prg_name = 'Atropos'
    file_type_name = 'Adapter trimmed files'
    ainfo = _per_sample_ainfo(
        out_dir, samples, suffixes, prg_name, file_type_name, bool(rs))

    return True, ainfo, ""
