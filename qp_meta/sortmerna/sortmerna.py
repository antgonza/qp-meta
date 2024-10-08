# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------


from os.path import join
from os import environ
from qp_meta.utils import (
    make_read_pairs_per_sample,
    _per_sample_ainfo, _generate_qiime_mapping_file)

DIR = environ["QC_SORTMERNA_DB_DP"]

RNA_REF_DB = ('--ref {0}smr_v4.3_default_db.fasta').format(DIR)


# resources per job
PPN = 10
MEMORY = '40g'
WALLTIME = '30:00:00'
FINISH_MEMORY = '48g'
FINISH_WALLTIME = '10:00:00'
MAX_RUNNING = 8


def generate_sortmerna_commands(forward_seqs, reverse_seqs, map_file,
                                out_dir, parameters):
    """Generates the Sortmerna commands

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
        The Sortmerna commands
    samples: list of tup
        list of 4-tuples with run prefix, sample name, fwd read fp, rev read fp

    Notes
    -----
    Currently this is requiring matched pairs in the make_read_pairs_per_sample
    step but implicitly allowing empty reverse reads in the actual command
    generation. This behavior may allow support of situations with empty
    reverse reads in some samples, for example after trimming and QC.
    """
    # matching filenames, samples, and run prefixes
    samples = make_read_pairs_per_sample(forward_seqs, reverse_seqs, map_file)

    cmds = []
    # --aligned {smr_r_op} --other {smr_nr_op}
    if reverse_seqs:
        template = (
            "sortmerna {ref_db} --reads {fwd} --reads {rev} --workdir {wkdir} "
            "--other --aligned --fastx --blast 1 --num_alignments 1 "
            "--threads {thrds} --paired_in --out2 -m 3988 --log")
    else:
        template = (
            "sortmerna {ref_db} --reads {fwd} --workdir {wkdir} "
            "--other --aligned --fastx --blast 1 --num_alignments 1 "
            "--threads {thrds} --out2 -m 3988 --log")

    arguments = {'thrds': PPN, 'ref_db': RNA_REF_DB}
    for run_prefix, sample, f_fp, r_fp in samples:
        arguments['wkdir'] = join(out_dir, run_prefix)
        arguments['fwd'] = f_fp
        arguments['rev'] = r_fp

        # arguments['smr_r_op_gz'] = arguments['smr_r_op'] + '.fastq.gz'
        # arguments['smr_nr_op_gz'] = arguments['smr_nr_op'] + '.fastq.gz'

        cmds.append(template.format(**arguments))

    return cmds, samples
    # In this version I have not added a summary file or sam file


def sortmerna(qclient, job_id, parameters, out_dir):
    """Run Sortmerna with the given parameters

    Parameters
    ----------
    qclient : tgp.qiita_client.QiitaClient
        The Qiita server client
    job_id : str
        The job id
    parameters : dict
        The parameter values
    out_dir : str
        The path to the job's output directory

    Returns
    -------
    bool, list, str
        The results of the job
    """
    qclient.update_job_step(
        job_id, "Step 3 of 4: Finishing SortMeRNA")

    ainfo = []

    # Generates 2 artifacts: one for the ribosomal
    # reads and other for the non-ribosomal reads
    samples = []
    with open(f'{out_dir}/{job_id}.samples.tsv') as f:
        for line in f.readlines():
            line = line.split()
            if len(line) == 3:
                line.append(None)
            samples.append(tuple(line))

    # checking the size of the first entry to see if there are reverse_seqs
    reverse = False
    if len(samples[0]) > 3:
        reverse = True

    suffixes = ['%s.nonribosomal.R1.fastq.gz', '%s.nonribosomal.R2.fastq.gz']
    file_type_name = 'Non-ribosomal reads'
    ainfo = [_per_sample_ainfo(
        out_dir, samples, suffixes, file_type_name, reverse)]

    suffixes = ['%s.ribosomal.R1.fastq.gz', '%s.ribosomal.R2.fastq.gz']
    file_type_name = 'Ribosomal reads'
    ainfo.append(_per_sample_ainfo(
        out_dir, samples, suffixes, file_type_name, reverse, add_logs=True))

    return True, ainfo, ""


def sortmerna_to_array(files, out_dir, params, prep_info, url, job_id):
    """Creates slurm files for submission of per sample bowtie2 and woltka
    """
    prep_file = _generate_qiime_mapping_file(prep_info, out_dir)

    reverse_seqs = []
    if 'raw_reverse_seqs' in files:
        reverse_seqs = files['raw_reverse_seqs']

    commands, samples = generate_sortmerna_commands(
        files['raw_forward_seqs'], reverse_seqs, prep_file, out_dir, params)

    # writing the job array details
    details_name = join(out_dir, 'sortmerna.array-details')
    with open(details_name, 'w') as details:
        details.write('\n'.join(commands))
    n_jobs = len(commands)

    # all the setup pieces
    lines = ['#!/bin/bash',
             '#SBATCH -p qiita',
             '#SBATCH --mail-user qiita.help@gmail.com',
             f'#SBATCH --job-name {job_id}',
             '#SBATCH -N 1',
             f'#SBATCH -n {PPN}',
             f'#SBATCH --time {WALLTIME}',
             f'#SBATCH --mem {MEMORY}',
             f'#SBATCH --output {out_dir}/{job_id}_%a.log',
             f'#SBATCH --error {out_dir}/{job_id}_%a.err',
             f'#SBATCH --array 1-{n_jobs}%{MAX_RUNNING}',
             'set -e',
             f'cd {out_dir}',
             f'{params["environment"]}',
             'date',  # start time
             'hostname',  # executing system
             'echo ${SLURM_JOBID} ${SLURM_ARRAY_TASK_ID}',
             'offset=${SLURM_ARRAY_TASK_ID}',
             'step=$(( $offset - 0 ))',
             f'cmd=$(head -n $step {details_name} | tail -n 1)',
             'eval $cmd',
             'set +e',
             'date']
    main_fp = join(out_dir, f'{job_id}.slurm')
    with open(main_fp, 'w') as job:
        job.write('\n'.join(lines))
        job.write('\n')

    # finish job
    lines = ['#!/bin/bash',
             '#SBATCH -p qiita',
             '#SBATCH --mail-user qiita.help@gmail.com',
             f'#SBATCH --job-name finish-{job_id}',
             '#SBATCH -N 1',
             '#SBATCH -n 1',
             f'#SBATCH --time {FINISH_WALLTIME}',
             f'#SBATCH --mem {FINISH_MEMORY}',
             f'#SBATCH --output {out_dir}/finish-{job_id}.log',
             f'#SBATCH --error {out_dir}/finish-{job_id}.err',
             'set -e',
             f'cd {out_dir}',
             f'{params["environment"]}',
             'date',  # start time
             'hostname',  # executing system
             'echo $SLURM_JOBID',
             f'finish_sortmerna {url} {job_id} {out_dir}\n'
             "date"]
    finish_fp = join(out_dir, f'{job_id}.finish.slurm')
    with open(finish_fp, 'w') as out:
        out.write('\n'.join(lines))
        out.write('\n')

    samples_fp = join(out_dir, f'{job_id}.samples.tsv')
    with open(samples_fp, 'w') as out:
        out.write('\n'.join(
            ['\t'.join([ss for ss in s if ss is not None]) for s in samples]))

    return main_fp, finish_fp, samples_fp
