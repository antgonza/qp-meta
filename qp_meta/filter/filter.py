# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from os.path import join
from tempfile import TemporaryDirectory
from qp_meta.utils import (
    _format_params, make_read_pairs_per_sample,
    _run_commands, _per_sample_ainfo)

BOWTIE2_PARAMS = {
    'x': 'Bowtie2 database to filter',
    'p': 'Number of threads'}


def generate_filter_commands(forward_seqs, reverse_seqs, map_file,
                             out_dir, temp_dir, parameters):
    """Generates the QC_Filter commands

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
        The QC_Filter commands
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

    param_string = _format_params(parameters, BOWTIE2_PARAMS)
    threads = parameters['Number of threads']

    for run_prefix, sample, f_fp, r_fp in samples:
        cmds.append('bowtie2 {params} --very-sensitive -1 {fwd_ip} -2 {rev_ip}'
                    ' | samtools view -f 12 -F 256 -b -o {bow_op}; '

                    'samtools sort -T {sample_path} -@ {thrds} -n -o {sam_op} '
                    '{sam_un_op}; '

                    'bedtools bamtofastq -i {sam_op} -fq {bedtools_op_one} '
                    '-fq2 {bedtools_op_two}; '
                    'pigz -p {thrds} -c {bedtools_op_one} > {gz_op_one}; '
                    'pigz -p {thrds} -c {bedtools_op_two} > {gz_op_two};'

                    .format(params=param_string, thrds=threads,
                            fwd_ip=f_fp, rev_ip=r_fp,
                            bow_op=join(temp_dir, '%s.unsorted.bam'
                                        % run_prefix),
                            sample_path=join(temp_dir, '%s' % run_prefix),
                            sam_op=join(temp_dir, '%s.bam' % run_prefix),
                            sam_un_op=join(temp_dir,
                                           '%s.unsorted.bam' % run_prefix),
                            bedtools_op_one=join(temp_dir,
                                                 '%s.R1.fastq'
                                                 % run_prefix),
                            bedtools_op_two=join(temp_dir,
                                                 '%s.R2.fastq'
                                                 % run_prefix),
                            gz_op_one=join(out_dir,
                                           '%s.R1.fastq.gz'
                                           % run_prefix),
                            gz_op_two=join(out_dir,
                                           '%s.R2.fastq.gz'
                                           % run_prefix)))

    return cmds, samples


def filter(qclient, job_id, parameters, out_dir):
    """Run filtering using Bowtie2 with the given parameters

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
    # Step 1 get the rest of the information need to run Bowtie2
    qclient.update_job_step(job_id, "Step 1 of 4: Collecting information")
    artifact_id = parameters['input']
    del parameters['input']

    # Get the artifact filepath information
    artifact_info = qclient.get("/qiita_db/artifacts/%s/" % artifact_id)
    fps = artifact_info['files']

    # Get the artifact metadata
    prep_info = qclient.get('/qiita_db/prep_template/%s/'
                            % artifact_info['prep_information'][0])
    prep_file = prep_info['prep-file']

    # Step 2 generating command
    qclient.update_job_step(job_id, "Step 2 of 4: Generating"
                                    " QC_Filter commands")
    # Creating temporary directory for intermediate files
    with TemporaryDirectory(dir=out_dir, prefix='filter_') as temp_dir:
        rs = fps['raw_reverse_seqs'] if 'raw_reverse_seqs' in fps else []
        commands, samples = generate_filter_commands(fps['raw_forward_seqs'],
                                                     rs, prep_file, out_dir,
                                                     temp_dir, parameters)

        # Step 3 execute filtering command
        len_cmd = len(commands)
        msg = "Step 3 of 4: Executing QC_Filter job (%d/{0})".format(len_cmd)
        success, msg = _run_commands(
            qclient, job_id, commands, msg, 'QC_Filter')
        if not success:
            return False, None, msg

    # Step 4 generating artifacts
    msg = "Step 4 of 4: Generating new artifacts (%d/{0})".format(len_cmd)
    suffixes = ['%s.R1.fastq.gz',
                '%s.R2.fastq.gz']
    prg_name = 'Filtering'
    file_type_name = 'Filtered files'
    ainfo = _per_sample_ainfo(
        out_dir, samples, suffixes, prg_name, file_type_name, bool(rs))

    return True, ainfo, ""
