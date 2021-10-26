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
    'a': 'Number of threads',
    'm': 'Memory'}


class QC_SortmernaTests(PluginTestCase):
    maxDiff = None

    def setUp(self):
        plugin("https://localhost:21174", 'register', 'ignored')

        self.params = {
                       'Output blast format': '1',
                       'Number of alignments': '1',
                       'Memory': '3988',
                       'Number of threads': '5'
        }
        self._clean_up_files = []

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def test_format_sortmerna_params(self):
        obs = _format_params(self.params, SORTMERNA_PARAMS)
        exp = (
               '-a 5 '
               '--blast 1 '
               '-m 3988 '
               '--num_alignments 1'
               )
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
                'command': dumps(['qp-meta', '2021.11', 'Sortmerna v2.1b']),
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
        artifact_info = self.qclient.get("/qiita_db/artifacts/%s/" % aid)

        # Get the artifact metadata
        prep_info = self.qclient.get('/qiita_db/prep_template/%s/' % pid)
        prep_file = prep_info['prep-file']

        # extra parameters
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)
        url = 'this-is-my-url'
        sdb = environ['QC_SORTMERNA_DB_DP']

        main_qsub_fp, finish_qsub_fp, samples_fp = sortmerna_to_array(
            artifact_info['files'], out_dir, self.params, prep_file,
            url, job_id)

        od = partial(join, out_dir)
        self.assertEqual(od(f'{job_id}.qsub'), main_qsub_fp)
        self.assertEqual(od(f'{job_id}.finish.qsub'), finish_qsub_fp)
        self.assertEqual(od(f'{job_id}.samples.tsv'), samples_fp)

        with open(main_qsub_fp) as f:
            main_qsub = f.readlines()
        with open(finish_qsub_fp) as f:
            finish_qsub = f.readlines()
        with open(samples_fp) as f:
            samples = f.readlines()

        exp_main_qsub = [
            '#!/bin/bash\n',
            '#PBS -M qiita.help@gmail.com\n',
            f'#PBS -N {job_id}\n',
            '#PBS -l nodes=1:ppn=10\n',
            '#PBS -l walltime=30:00:00\n',
            '#PBS -l mem=40g\n',
            f'#PBS -o {out_dir}/{job_id}_${{PBS_ARRAYID}}.log\n',
            f'#PBS -e {out_dir}/{job_id}_${{PBS_ARRAYID}}.err\n',
            '#PBS -t 1-4%8\n',
            '#PBS -l epilogue=/home/qiita/qiita-epilogue.sh\n',
            'set -e\n',
            f'cd {out_dir}\n',
            f'export QC_SORTMERNA_DB_DP={sdb}; '
            'source ~/.bash_profile; conda activate qp-meta\n',
            'date\n',
            'hostname\n',
            'echo ${PBS_JOBID} ${PBS_ARRAYID}\n',
            'offset=${PBS_ARRAYID}\n',
            'step=$(( $offset - 0 ))\n',
            f'cmd=$(head -n $step {out_dir}/sortmerna.array-details | '
            'tail -n 1)\n',
            'eval $cmd\n',
            'set +e\n',
            'date\n']
        self.assertEqual(main_qsub, exp_main_qsub)

        exp_finish_qsub = [
            '#!/bin/bash\n',
            '#PBS -M qiita.help@gmail.com\n',
            f'#PBS -N finish-{job_id}\n',
            '#PBS -l nodes=1:ppn=1\n',
            '#PBS -l walltime=10:00:00\n',
            '#PBS -l mem=48g\n',
            f'#PBS -o {out_dir}/finish-{job_id}.log\n',
            f'#PBS -e {out_dir}/finish-{job_id}.err\n',
            '#PBS -l epilogue=/home/qiita/qiita-epilogue.sh\n',
            'set -e\n',
            f'cd {out_dir}\n',
            f'export QC_SORTMERNA_DB_DP={sdb}; '
            'source ~/.bash_profile; conda activate qp-meta\n',
            'date\n',
            'hostname\n',
            'echo $PBS_JOBID\n',
            f'finish_sortmerna {url} {job_id} {out_dir}\n',
            'date\n']
        self.assertEqual(finish_qsub, exp_finish_qsub)

        adir = dirname(artifact_info['files']['raw_forward_seqs'][0])
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
            f'unpigz -p 5 -c {adir}/S22205_S104_L001_R1_001.fastq.gz > '
            f'{out_dir}/S22205_S104_L001_R1_001.fastq && sortmerna --ref '
            f'{RNA_REF_DB} --reads {out_dir}/S22205_S104_L001_R1_001.fastq '
            f'--aligned {out_dir}/S22205_S104.ribosomal.R1 --other '
            f'{out_dir}/S22205_S104.nonribosomal.R1 --fastx -a 5 --blast 1 '
            '-m 3988 --num_alignments 1 && pigz -p 5 -c '
            f'{out_dir}/S22205_S104.ribosomal.R1.fastq > '
            f'{out_dir}/S22205_S104.ribosomal.R1.fastq.gz && pigz -p 5 -c '
            f'{out_dir}/S22205_S104.nonribosomal.R1.fastq > '
            f'{out_dir}/S22205_S104.nonribosomal.R1.fastq.gz;\n',
            f'unpigz -p 5 -c {adir}/S22205_S104_L001_R2_001.fastq.gz > '
            f'{out_dir}/S22205_S104_L001_R2_001.fastq && sortmerna --ref '
            f'{RNA_REF_DB} --reads {out_dir}/S22205_S104_L001_R2_001.fastq '
            f'--aligned {out_dir}/S22205_S104.ribosomal.R2 --other '
            f'{out_dir}/S22205_S104.nonribosomal.R2 --fastx -a 5 --blast 1 '
            '-m 3988 --num_alignments 1 && pigz -p 5 -c '
            f'{out_dir}/S22205_S104.ribosomal.R2.fastq > '
            f'{out_dir}/S22205_S104.ribosomal.R2.fastq.gz && pigz -p 5 -c '
            f'{out_dir}/S22205_S104.nonribosomal.R2.fastq > '
            f'{out_dir}/S22205_S104.nonribosomal.R2.fastq.gz;\n',
            f'unpigz -p 5 -c {adir}/S22282_S102_L001_R1_001.fastq.gz > '
            f'{out_dir}/S22282_S102_L001_R1_001.fastq && sortmerna --ref '
            f'{RNA_REF_DB} --reads {out_dir}/S22282_S102_L001_R1_001.fastq '
            f'--aligned {out_dir}/S22282_S102.ribosomal.R1 --other '
            f'{out_dir}/S22282_S102.nonribosomal.R1 --fastx -a 5 --blast 1 '
            '-m 3988 --num_alignments 1 && pigz -p 5 -c '
            f'{out_dir}/S22282_S102.ribosomal.R1.fastq > '
            f'{out_dir}/S22282_S102.ribosomal.R1.fastq.gz && pigz -p 5 -c '
            f'{out_dir}/S22282_S102.nonribosomal.R1.fastq > '
            f'{out_dir}/S22282_S102.nonribosomal.R1.fastq.gz;\n',
            f'unpigz -p 5 -c {adir}/S22282_S102_L001_R2_001.fastq.gz > '
            f'{out_dir}/S22282_S102_L001_R2_001.fastq && sortmerna --ref '
            f'{RNA_REF_DB} --reads {out_dir}/S22282_S102_L001_R2_001.fastq '
            f'--aligned {out_dir}/S22282_S102.ribosomal.R2 --other '
            f'{out_dir}/S22282_S102.nonribosomal.R2 --fastx -a 5 --blast 1 '
            '-m 3988 --num_alignments 1 && pigz -p 5 -c '
            f'{out_dir}/S22282_S102.ribosomal.R2.fastq > '
            f'{out_dir}/S22282_S102.ribosomal.R2.fastq.gz && pigz -p 5 -c '
            f'{out_dir}/S22282_S102.nonribosomal.R2.fastq > '
            f'{out_dir}/S22282_S102.nonribosomal.R2.fastq.gz;']
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
                 'raw_reverse_seqs')]
        nribo = [(f'{out_dir}/S22282_S102.nonribosomal.R1.fastq.gz',
                  'raw_forward_seqs'),
                 (f'{out_dir}/S22282_S102.nonribosomal.R2.fastq.gz',
                  'raw_reverse_seqs'),
                 (f'{out_dir}/S22205_S104.nonribosomal.R1.fastq.gz',
                  'raw_forward_seqs'),
                 (f'{out_dir}/S22205_S104.nonribosomal.R2.fastq.gz',
                  'raw_reverse_seqs')]
        for f in ribo + nribo:
            copyfile(infile, f[0])
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
        artifact_info = self.qclient.get("/qiita_db/artifacts/%s/" % aid)

        # Get the artifact metadata
        prep_info = self.qclient.get('/qiita_db/prep_template/%s/' % pid)
        prep_file = prep_info['prep-file']

        # extra parameters
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)
        url = 'this-is-my-url'
        sdb = environ['QC_SORTMERNA_DB_DP']

        main_qsub_fp, finish_qsub_fp, samples_fp = sortmerna_to_array(
            artifact_info['files'], out_dir, self.params, prep_file,
            url, job_id)

        od = partial(join, out_dir)
        self.assertEqual(od(f'{job_id}.qsub'), main_qsub_fp)
        self.assertEqual(od(f'{job_id}.finish.qsub'), finish_qsub_fp)
        self.assertEqual(od(f'{job_id}.samples.tsv'), samples_fp)

        with open(main_qsub_fp) as f:
            main_qsub = f.readlines()
        with open(finish_qsub_fp) as f:
            finish_qsub = f.readlines()
        with open(samples_fp) as f:
            samples = f.readlines()

        exp_main_qsub = [
            '#!/bin/bash\n',
            '#PBS -M qiita.help@gmail.com\n',
            f'#PBS -N {job_id}\n',
            '#PBS -l nodes=1:ppn=10\n',
            '#PBS -l walltime=30:00:00\n',
            '#PBS -l mem=40g\n',
            f'#PBS -o {out_dir}/{job_id}_${{PBS_ARRAYID}}.log\n',
            f'#PBS -e {out_dir}/{job_id}_${{PBS_ARRAYID}}.err\n',
            '#PBS -t 1-2%8\n',
            '#PBS -l epilogue=/home/qiita/qiita-epilogue.sh\n',
            'set -e\n',
            f'cd {out_dir}\n',
            f'export QC_SORTMERNA_DB_DP={sdb}; '
            'source ~/.bash_profile; conda activate qp-meta\n',
            'date\n',
            'hostname\n',
            'echo ${PBS_JOBID} ${PBS_ARRAYID}\n',
            'offset=${PBS_ARRAYID}\n',
            'step=$(( $offset - 0 ))\n',
            f'cmd=$(head -n $step {out_dir}/sortmerna.array-details | '
            'tail -n 1)\n',
            'eval $cmd\n',
            'set +e\n',
            'date\n']
        self.assertEqual(main_qsub, exp_main_qsub)

        exp_finish_qsub = [
            '#!/bin/bash\n',
            '#PBS -M qiita.help@gmail.com\n',
            f'#PBS -N finish-{job_id}\n',
            '#PBS -l nodes=1:ppn=1\n',
            '#PBS -l walltime=10:00:00\n',
            '#PBS -l mem=48g\n',
            f'#PBS -o {out_dir}/finish-{job_id}.log\n',
            f'#PBS -e {out_dir}/finish-{job_id}.err\n',
            '#PBS -l epilogue=/home/qiita/qiita-epilogue.sh\n',
            'set -e\n',
            f'cd {out_dir}\n',
            f'export QC_SORTMERNA_DB_DP={sdb}; '
            'source ~/.bash_profile; conda activate qp-meta\n',
            'date\n',
            'hostname\n',
            'echo $PBS_JOBID\n',
            f'finish_sortmerna {url} {job_id} {out_dir}\n',
            'date\n']
        self.assertEqual(finish_qsub, exp_finish_qsub)

        adir = dirname(artifact_info['files']['raw_forward_seqs'][0])
        exp_samples = [
            'S22205_S104\t1.SKB8.640193\t'
            f'{adir}/S22205_S104_L001_R1_001.fastq.gz\n',
            'S22282_S102\t1.SKD8.640184\t'
            f'{adir}/S22282_S102_L001_R1_001.fastq.gz']
        self.assertEqual(exp_samples, samples)

        with open(f'{out_dir}/sortmerna.array-details') as f:
            details = f.readlines()
        exp_details = [
            f'unpigz -p 5 -c {adir}/S22205_S104_L001_R1_001.fastq.gz > '
            f'{out_dir}/S22205_S104_L001_R1_001.fastq && sortmerna --ref '
            f'{RNA_REF_DB} --reads {out_dir}/S22205_S104_L001_R1_001.fastq '
            f'--aligned {out_dir}/S22205_S104.ribosomal.R1 --other '
            f'{out_dir}/S22205_S104.nonribosomal.R1 --fastx -a 5 --blast 1 '
            '-m 3988 --num_alignments 1 && pigz -p 5 -c '
            f'{out_dir}/S22205_S104.ribosomal.R1.fastq > '
            f'{out_dir}/S22205_S104.ribosomal.R1.fastq.gz && pigz -p 5 -c '
            f'{out_dir}/S22205_S104.nonribosomal.R1.fastq > '
            f'{out_dir}/S22205_S104.nonribosomal.R1.fastq.gz;\n',
            f'unpigz -p 5 -c {adir}/S22282_S102_L001_R1_001.fastq.gz > '
            f'{out_dir}/S22282_S102_L001_R1_001.fastq && sortmerna --ref '
            f'{RNA_REF_DB} --reads {out_dir}/S22282_S102_L001_R1_001.fastq '
            f'--aligned {out_dir}/S22282_S102.ribosomal.R1 --other '
            f'{out_dir}/S22282_S102.nonribosomal.R1 --fastx -a 5 --blast 1 '
            '-m 3988 --num_alignments 1 && pigz -p 5 -c '
            f'{out_dir}/S22282_S102.ribosomal.R1.fastq > '
            f'{out_dir}/S22282_S102.ribosomal.R1.fastq.gz && pigz -p 5 -c '
            f'{out_dir}/S22282_S102.nonribosomal.R1.fastq > '
            f'{out_dir}/S22282_S102.nonribosomal.R1.fastq.gz;']
        self.assertEqual(exp_details, details)

        # making sure it finishes correctly
        infile = f'{adir}/S22205_S104_L001_R1_001.fastq.gz'
        ribo = [(f'{out_dir}/S22282_S102.ribosomal.R1.fastq.gz',
                 'raw_forward_seqs'),
                (f'{out_dir}/S22205_S104.ribosomal.R1.fastq.gz',
                 'raw_forward_seqs')]
        nribo = [(f'{out_dir}/S22282_S102.nonribosomal.R1.fastq.gz',
                  'raw_forward_seqs'),
                 (f'{out_dir}/S22205_S104.nonribosomal.R1.fastq.gz',
                  'raw_forward_seqs')]
        for f in ribo + nribo:
            copyfile(infile, f[0])
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
