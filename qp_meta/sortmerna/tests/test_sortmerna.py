# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import main
from os import remove, makedirs, environ
from os.path import exists, isdir, join, dirname
from shutil import rmtree, copyfile
from tempfile import mkdtemp
from json import dumps
from functools import partial

from qiita_client.testing import PluginTestCase
from qiita_client import ArtifactInfo

from qp_meta import plugin
from qp_meta.sortmerna.sortmerna import (
    sortmerna_to_array, sortmerna, RNA_REF_DB)
from qp_meta.utils import (
    _format_params, _per_sample_ainfo)

SORTMERNA_PARAMS = {
    'blast': 'Output blast format',
    'num_alignments': 'Number of alignments',
    'm': 'Memory'}


class QC_SortmernaTests(PluginTestCase):
    maxDiff = None

    def setUp(self):
        plugin("https://localhost:21174", 'register', 'ignored')

        self.params = {}
        self._clean_up_files = []

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def test_format_sortmerna_params(self):
        obs = _format_params(
            {'Output blast format': '1',
             'Number of alignments': '1',
             'Memory': '3988'}, SORTMERNA_PARAMS)
        exp = ('--blast 1 '
               '-m 3988 '
               '--num_alignments 1')
        self.assertEqual(obs, exp)

    def _helper_tester(self, prep_info_dict, just_forward=True):
        data = {'prep_info': dumps(prep_info_dict),
                # magic #1 = testing study
                'study': 1,
                'data_type': 'Metagenomic'}
        pid = self.qclient.post('/apitest/prep_template/', data=data)['prep']

        # inserting artifacts
        in_dir = mkdtemp()
        self._clean_up_files.append(in_dir)

        self.source_dir = 'qp_meta/sortmerna/tests/support_files/'
        fp1_1 = join(in_dir, 'S22205_S104_L001_R1_001.fastq.gz')
        fp2_1 = join(in_dir, 'S22282_S102_L001_R1_001.fastq.gz')
        copyfile(f'{self.source_dir}/S22205_S104_L001_R1_001.fastq.gz', fp1_1)
        copyfile(f'{self.source_dir}/S22282_S102_L001_R1_001.fastq.gz', fp2_1)
        filepaths = [(fp1_1, 'raw_forward_seqs'), (fp2_1, 'raw_forward_seqs')]
        if not just_forward:
            fp1_2 = join(in_dir, 'S22205_S104_L001_R2_001.fastq.gz')
            fp2_2 = join(in_dir, 'S22282_S102_L001_R2_001.fastq.gz')
            copyfile(
                f'{self.source_dir}/S22205_S104_L001_R2_001.fastq.gz', fp1_2)
            copyfile(
                f'{self.source_dir}/S22282_S102_L001_R2_001.fastq.gz', fp2_2)
            filepaths.extend(
                [(fp1_2, 'raw_reverse_seqs'), (fp2_2, 'raw_reverse_seqs')])

        data = {
            'filepaths': dumps(filepaths),
            'type': "per_sample_FASTQ",
            'name': "Test artifact",
            'prep': pid}
        aid = self.qclient.post('/apitest/artifact/', data=data)['artifact']

        self.params['input'] = aid

        data = {'user': 'demo@microbio.me',
                'command': dumps(['qp-meta', '2024.10', 'SortMeRNA v4.3.7']),
                'status': 'running',
                'parameters': dumps(self.params)}
        job_id = self.qclient.post(
            '/apitest/processing_job/', data=data)['job']

        return pid, aid, job_id

    def test_generate_sortmerna_analysis_commands_forward_reverse(self):
        prep_info_dict = {
            'SKB8.640193': {'run_prefix': 'S22205_S104'},
            'SKD8.640184': {'run_prefix': 'S22282_S102'}}
        pid, aid, job_id = self._helper_tester(prep_info_dict, False)

        # adding extra parameters
        self.params['environment'] = environ["ENVIRONMENT"]

        # Get the artifact filepath information
        files, prep = self.qclient.artifact_and_preparation_files(aid)
        fps = {'raw_forward_seqs': [], 'raw_reverse_seqs': []}
        for sn, fs in files.items():
            fps['raw_forward_seqs'].append(fs[0]['filepath'])
            if fs[1]:
                fps['raw_reverse_seqs'].append(fs[1]['filepath'])

        # extra parameters
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)
        url = 'this-is-my-url'

        main_fp, finish_fp, samples_fp = sortmerna_to_array(
            fps, out_dir, self.params, prep, url, job_id)

        od = partial(join, out_dir)
        self.assertEqual(od(f'{job_id}.slurm'), main_fp)
        self.assertEqual(od(f'{job_id}.finish.slurm'), finish_fp)
        self.assertEqual(od(f'{job_id}.samples.tsv'), samples_fp)

        with open(main_fp) as f:
            main = f.readlines()
        with open(finish_fp) as f:
            finish = f.readlines()
        with open(samples_fp) as f:
            samples = f.readlines()

        exp_main = [
            '#!/bin/bash\n',
            '#SBATCH -p qiita\n',
            '#SBATCH --mail-user qiita.help@gmail.com\n',
            f'#SBATCH --job-name {job_id}\n',
            '#SBATCH -N 1\n',
            '#SBATCH -n 10\n',
            '#SBATCH --time 30:00:00\n',
            '#SBATCH --mem 40g\n',
            f'#SBATCH --output {out_dir}/{job_id}_%a.log\n',
            f'#SBATCH --error {out_dir}/{job_id}_%a.err\n',
            '#SBATCH --array 1-2%8\n',
            'set -e\n',
            f'cd {out_dir}\n',
            f'{self.params["environment"]}\n',
            'date\n',
            'hostname\n',
            'echo ${SLURM_JOBID} ${SLURM_ARRAY_TASK_ID}\n',
            'offset=${SLURM_ARRAY_TASK_ID}\n',
            'step=$(( $offset - 0 ))\n',
            f'cmd=$(head -n $step {out_dir}/sortmerna.array-details | '
            'tail -n 1)\n',
            'eval $cmd\n',
            'set +e\n',
            'date\n']
        self.assertEqual(main, exp_main)

        exp_finish = [
            '#!/bin/bash\n',
            '#SBATCH -p qiita\n',
            '#SBATCH --mail-user qiita.help@gmail.com\n',
            f'#SBATCH --job-name finish-{job_id}\n',
            '#SBATCH -N 1\n',
            '#SBATCH -n 1\n',
            '#SBATCH --time 10:00:00\n',
            '#SBATCH --mem 1g\n',
            f'#SBATCH --output {out_dir}/finish-{job_id}.log\n',
            f'#SBATCH --error {out_dir}/finish-{job_id}.err\n',
            'set -e\n',
            f'cd {out_dir}\n',
            f'{self.params["environment"]}\n',
            'date\n',
            'hostname\n',
            'echo $SLURM_JOBID\n',
            f'finish_sortmerna {url} {job_id} {out_dir}\n',
            'date\n']
        self.assertEqual(finish, exp_finish)

        adir = dirname(fps['raw_forward_seqs'][0])
        exp_samples = [
            'S22205_S104\t1.SKB8.640193\t'
            f'{adir}/S22205_S104_L001_R1_001.fastq.gz\t'
            f'{adir}/S22205_S104_L001_R2_001.fastq.gz\n',
            'S22282_S102\t1.SKD8.640184\t'
            f'{adir}/S22282_S102_L001_R1_001.fastq.gz\t'
            f'{adir}/S22282_S102_L001_R2_001.fastq.gz']
        self.assertEqual(exp_samples, samples)

        with open(f'{out_dir}/sortmerna.array-details') as f:
            details = f.readlines()
        exp_details = [
            f'sortmerna {RNA_REF_DB} '
            f'--reads {adir}/S22205_S104_L001_R1_001.fastq.gz '
            f'--reads {adir}/S22205_S104_L001_R2_001.fastq.gz '
            f'--workdir {out_dir}/S22205_S104 --other --aligned --fastx '
            f'--blast 1 --num_alignments 1 --threads 10 --paired_in --out2 '
            '-index 0; '
            f"mv {out_dir}/S22205_S104/out/aligned.log "
            f"{out_dir}/S22205_S104.log; "
            f'mv {out_dir}/S22205_S104/out/aligned_fwd.fq.gz '
            f'{out_dir}/S22205_S104.ribosomal.R1.fastq.gz; '
            f'mv {out_dir}/S22205_S104/out/aligned_rev.fq.gz '
            f'{out_dir}/S22205_S104.ribosomal.R2.fastq.gz; '
            f'mv {out_dir}/S22205_S104/out/other_fwd.fq.gz '
            f'{out_dir}/S22205_S104.nonribosomal.R1.fastq.gz; '
            f'mv {out_dir}/S22205_S104/out/other_rev.fq.gz '
            f'{out_dir}/S22205_S104.nonribosomal.R2.fastq.gz; \n',
            f'sortmerna {RNA_REF_DB} '
            f'--reads {adir}/S22282_S102_L001_R1_001.fastq.gz '
            f'--reads {adir}/S22282_S102_L001_R2_001.fastq.gz '
            f'--workdir {out_dir}/S22282_S102 --other --aligned --fastx '
            f'--blast 1 --num_alignments 1 --threads 10 --paired_in --out2 '
            '-index 0; '
            f"mv {out_dir}/S22282_S102/out/aligned.log "
            f"{out_dir}/S22282_S102.log; "
            f'mv {out_dir}/S22282_S102/out/aligned_fwd.fq.gz '
            f'{out_dir}/S22282_S102.ribosomal.R1.fastq.gz; '
            f'mv {out_dir}/S22282_S102/out/aligned_rev.fq.gz '
            f'{out_dir}/S22282_S102.ribosomal.R2.fastq.gz; '
            f'mv {out_dir}/S22282_S102/out/other_fwd.fq.gz '
            f'{out_dir}/S22282_S102.nonribosomal.R1.fastq.gz; '
            f'mv {out_dir}/S22282_S102/out/other_rev.fq.gz '
            f'{out_dir}/S22282_S102.nonribosomal.R2.fastq.gz; ']
        self.assertEqual(exp_details, details)

        # making sure it finishes correctly
        infile = f'{adir}/S22205_S104_L001_R1_001.fastq.gz'
        ribo = [(f'{out_dir}/S22282_S102.ribosomal.R1.fastq.gz',
                 'raw_forward_seqs'),
                (f'{out_dir}/S22282_S102.ribosomal.R2.fastq.gz',
                 'raw_reverse_seqs'),
                (f'{out_dir}/S22205_S104.ribosomal.R1.fastq.gz',
                 'raw_forward_seqs'),
                (f'{out_dir}/S22205_S104.ribosomal.R2.fastq.gz',
                 'raw_reverse_seqs'),
                (f'{out_dir}/artifact.logs.tgz', 'log')]
        nribo = [(f'{out_dir}/S22282_S102.nonribosomal.R1.fastq.gz',
                  'raw_forward_seqs'),
                 (f'{out_dir}/S22282_S102.nonribosomal.R2.fastq.gz',
                  'raw_reverse_seqs'),
                 (f'{out_dir}/S22205_S104.nonribosomal.R1.fastq.gz',
                  'raw_forward_seqs'),
                 (f'{out_dir}/S22205_S104.nonribosomal.R2.fastq.gz',
                  'raw_reverse_seqs')]
        for f in ribo + nribo:
            if f[1] == 'log':
                continue
            copyfile(infile, f[0])
            copyfile(infile, f[0].replace('.ribosomal.R1.fastq.gz', '.log'))
        success, ainfo, msg = sortmerna(
            self.qclient, job_id, self.params, out_dir)

        self.assertEqual("", msg)
        self.assertTrue(success)

        exp = [
            ArtifactInfo('Non-ribosomal reads', 'per_sample_FASTQ', nribo),
            ArtifactInfo('Ribosomal reads', 'per_sample_FASTQ', ribo)]
        self.assertCountEqual(ainfo, exp)

    def test_generate_sortmerna_analysis_commands_forward(self):
        prep_info_dict = {
            'SKB8.640193': {'run_prefix': 'S22205_S104'},
            'SKD8.640184': {'run_prefix': 'S22282_S102'}}
        pid, aid, job_id = self._helper_tester(prep_info_dict)

        # adding extra parameters
        self.params['environment'] = environ["ENVIRONMENT"]

        # Get the artifact filepath information
        files, prep = self.qclient.artifact_and_preparation_files(aid)
        fps = {'raw_forward_seqs': [], 'raw_reverse_seqs': []}
        for sn, fs in files.items():
            fps['raw_forward_seqs'].append(fs[0]['filepath'])
            if fs[1]:
                fps['raw_reverse_seqs'].append(fs[1]['filepath'])

        # extra parameters
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)
        url = 'this-is-my-url'

        main_fp, finish_fp, samples_fp = sortmerna_to_array(
            fps, out_dir, self.params, prep, url, job_id)

        od = partial(join, out_dir)
        self.assertEqual(od(f'{job_id}.slurm'), main_fp)
        self.assertEqual(od(f'{job_id}.finish.slurm'), finish_fp)
        self.assertEqual(od(f'{job_id}.samples.tsv'), samples_fp)

        with open(main_fp) as f:
            main = f.readlines()
        with open(finish_fp) as f:
            finish = f.readlines()
        with open(samples_fp) as f:
            samples = f.readlines()

        exp_main = [
            '#!/bin/bash\n',
            '#SBATCH -p qiita\n',
            '#SBATCH --mail-user qiita.help@gmail.com\n',
            f'#SBATCH --job-name {job_id}\n',
            '#SBATCH -N 1\n',
            '#SBATCH -n 10\n',
            '#SBATCH --time 30:00:00\n',
            '#SBATCH --mem 40g\n',
            f'#SBATCH --output {out_dir}/{job_id}_%a.log\n',
            f'#SBATCH --error {out_dir}/{job_id}_%a.err\n',
            '#SBATCH --array 1-2%8\n',
            'set -e\n',
            f'cd {out_dir}\n',
            f'{self.params["environment"]}\n',
            'date\n',
            'hostname\n',
            'echo ${SLURM_JOBID} ${SLURM_ARRAY_TASK_ID}\n',
            'offset=${SLURM_ARRAY_TASK_ID}\n',
            'step=$(( $offset - 0 ))\n',
            f'cmd=$(head -n $step {out_dir}/sortmerna.array-details | '
            'tail -n 1)\n',
            'eval $cmd\n',
            'set +e\n',
            'date\n']
        self.assertEqual(main, exp_main)

        exp_finish = [
            '#!/bin/bash\n',
            '#SBATCH -p qiita\n',
            '#SBATCH --mail-user qiita.help@gmail.com\n',
            f'#SBATCH --job-name finish-{job_id}\n',
            '#SBATCH -N 1\n',
            '#SBATCH -n 1\n',
            '#SBATCH --time 10:00:00\n',
            '#SBATCH --mem 1g\n',
            f'#SBATCH --output {out_dir}/finish-{job_id}.log\n',
            f'#SBATCH --error {out_dir}/finish-{job_id}.err\n',
            'set -e\n',
            f'cd {out_dir}\n',
            f'{self.params["environment"]}\n',
            'date\n',
            'hostname\n',
            'echo $SLURM_JOBID\n',
            f'finish_sortmerna {url} {job_id} {out_dir}\n',
            'date\n']
        self.assertEqual(finish, exp_finish)

        adir = dirname(fps['raw_forward_seqs'][0])
        exp_samples = [
            'S22205_S104\t1.SKB8.640193\t'
            f'{adir}/S22205_S104_L001_R1_001.fastq.gz\n',
            'S22282_S102\t1.SKD8.640184\t'
            f'{adir}/S22282_S102_L001_R1_001.fastq.gz']
        self.assertEqual(exp_samples, samples)

        with open(f'{out_dir}/sortmerna.array-details') as f:
            details = f.readlines()
        exp_details = [
            f'sortmerna {RNA_REF_DB} '
            f'--reads {adir}/S22205_S104_L001_R1_001.fastq.gz '
            f'--workdir {out_dir}/S22205_S104 --other --aligned --fastx '
            f'--blast 1 --num_alignments 1 --threads 10 '
            '--out2 -index 0; '
            f'mv {out_dir}/S22205_S104/out/aligned.log '
            f'{out_dir}/S22205_S104.log; '
            f'mv {out_dir}/S22205_S104/out/aligned_fwd.fq.gz '
            f'{out_dir}/S22205_S104.ribosomal.R1.fastq.gz; '
            f'mv {out_dir}/S22205_S104/out/other_fwd.fq.gz '
            f'{out_dir}/S22205_S104.nonribosomal.R1.fastq.gz; \n',
            f'sortmerna {RNA_REF_DB} '
            f'--reads {adir}/S22282_S102_L001_R1_001.fastq.gz '
            f'--workdir {out_dir}/S22282_S102 --other --aligned --fastx '
            f'--blast 1 --num_alignments 1 --threads 10 '
            '--out2 -index 0; '
            f"mv {out_dir}/S22282_S102/out/aligned.log "
            f"{out_dir}/S22282_S102.log; "
            f'mv {out_dir}/S22282_S102/out/aligned_fwd.fq.gz '
            f'{out_dir}/S22282_S102.ribosomal.R1.fastq.gz; '
            f'mv {out_dir}/S22282_S102/out/other_fwd.fq.gz '
            f'{out_dir}/S22282_S102.nonribosomal.R1.fastq.gz; ']
        self.assertEqual(exp_details, details)

        # making sure it finishes correctly
        infile = f'{adir}/S22205_S104_L001_R1_001.fastq.gz'
        ribo = [(f'{out_dir}/S22282_S102.ribosomal.R1.fastq.gz',
                 'raw_forward_seqs'),
                (f'{out_dir}/S22205_S104.ribosomal.R1.fastq.gz',
                 'raw_forward_seqs'),
                (f'{out_dir}/artifact.logs.tgz', 'log')]
        nribo = [(f'{out_dir}/S22282_S102.nonribosomal.R1.fastq.gz',
                  'raw_forward_seqs'),
                 (f'{out_dir}/S22205_S104.nonribosomal.R1.fastq.gz',
                  'raw_forward_seqs')]
        for f in ribo + nribo:
            if f[1] == 'log':
                continue
            copyfile(infile, f[0])
            copyfile(infile, f[0].replace('.ribosomal.R1.fastq.gz', '.log'))
        success, ainfo, msg = sortmerna(
            self.qclient, job_id, self.params, out_dir)

        self.assertEqual("", msg)
        self.assertTrue(success)

        exp = [
            ArtifactInfo('Non-ribosomal reads', 'per_sample_FASTQ', nribo),
            ArtifactInfo('Ribosomal reads', 'per_sample_FASTQ', ribo)]
        self.assertCountEqual(ainfo, exp)

    def test_per_sample_ainfo_error(self):
        in_dir = mkdtemp()
        self._clean_up_files.append(in_dir)
        makedirs(join(in_dir, 'sampleA'))
        makedirs(join(in_dir, 'sampleB'))

        # Paired-end
        with self.assertRaises(ValueError):
            _per_sample_ainfo(in_dir, (('sampleA', None, None, None),
                                       ('sampleB', None, None, None)), [],
                              'QC_Sortmerna Files', True)


if __name__ == '__main__':
    main()
