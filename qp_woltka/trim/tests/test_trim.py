# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import main
from os import close, remove, makedirs
from os.path import exists, isdir, join, dirname
from shutil import rmtree, copyfile
from tempfile import mkstemp, mkdtemp
from json import dumps
from functools import partial

from qiita_client.testing import PluginTestCase

from qp_woltka import plugin
from qp_woltka.trim.trim import (generate_trim_commands, trim)
from qp_woltka.utils import (
    _format_params, make_read_pairs_per_sample, _per_sample_ainfo)
import qp_woltka.trim as kd

ATROPOS_PARAMS = {
    'adapter': 'Fwd read adapter', 'A': 'Rev read adapter',
    'quality-cutoff': 'Trim low-quality bases',
    'minimum-length': 'Minimum trimmed read length',
    'pair-filter': 'Pair-end read required to match',
    'max-n': 'Maximum number of N bases in a read to keep it',
    'trim-n': 'Trim Ns on ends of reads', 'threads': 'Number of threads used',
    'nextseq-trim': 'NextSeq-specific quality trimming'}


class QC_TrimTests(PluginTestCase):
    maxDiff = None

    def setUp(self):
        plugin("https://localhost:21174", 'register', 'ignored')

        # to fully test the plugin we need to use the demo reference-db,
        # which is part of the regular trim install
        self.refdb_path = join(dirname(kd.__file__),
                               "tests/data/demo_bowtie2_db/demo_db.1.bt2")
        self.params = {
                       'Fwd read adapter': 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC',
                       'Rev read adapter': 'GATCGGAAGAGCGTCGTGTAGGGAAAGGAGTGT',
                       'Trim low-quality bases': '15',
                       'Minimum trimmed read length': '80',
                       'Pair-end read required to match': 'any',
                       'Maximum number of N bases in a read to keep it': '80',
                       'Trim Ns on ends of reads': True,
                       'NextSeq-specific quality trimming': False,
                       'Number of threads used': '5'
        }
        self._clean_up_files = []

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def test_format_trim_params(self):
        obs = _format_params(self.params, ATROPOS_PARAMS)
        exp = ('-A GATCGGAAGAGCGTCGTGTAGGGAAAGGAGTGT --adapter GATCGGAAGAGCACA'
               'CGTCTGAACTCCAGTCAC --max-n 80 --minimum-length 80 '
               '--pair-filter any --quality-cutoff 15 --threads 5 --trim-n')
        self.assertEqual(obs, exp)

    def test_make_read_pairs_per_sample_match_fwd_rev(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)

        fwd_fp = ['./folder/s3_S013_L001_R1.fastq.gz',
                  './folder/s2_S011_L001_R1.fastq.gz',
                  './folder/s1_S009_L001_R1.fastq.gz']

        rev_fp = ['./folder/s3_S013_L001_R2.fastq.gz',
                  './folder/s2_S011_L001_R2.fastq.gz',
                  './folder/s1_S009_L001_R2.fastq.gz']

        exp = [('s1', 'SKB8.640193', './folder/s1_S009_L001_R1.fastq.gz',
                './folder/s1_S009_L001_R2.fastq.gz'),
               ('s2', 'SKD8.640184', './folder/s2_S011_L001_R1.fastq.gz',
                './folder/s2_S011_L001_R2.fastq.gz'),
               ('s3', 'SKB7.640196', './folder/s3_S013_L001_R1.fastq.gz',
                './folder/s3_S013_L001_R2.fastq.gz')]

        obs = make_read_pairs_per_sample(fwd_fp, rev_fp, fp)

        self.assertEqual(obs, exp)

    def test_make_read_pairs_per_sample_match_fwd_only(self):
        # do we want to support forward-only yet?
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)

        fwd_fp = ['./folder/s3_S013_L001_R1.fastq.gz',
                  './folder/s2_S011_L001_R1.fastq.gz',
                  './folder/s1_S009_L001_R1.fastq.gz']

        rev_fp = []

        exp = [('s1', 'SKB8.640193', './folder/s1_S009_L001_R1.fastq.gz',
                None),
               ('s2', 'SKD8.640184', './folder/s2_S011_L001_R1.fastq.gz',
                None),
               ('s3', 'SKB7.640196', './folder/s3_S013_L001_R1.fastq.gz',
                None)]

        obs = make_read_pairs_per_sample(fwd_fp, rev_fp, fp)

        self.assertEqual(obs, exp)

    def test_make_read_pairs_per_sample_match_fwd_rev_diffnum(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)

        fwd_fp = ['./folder/s3_S013_L001_R1.fastq.gz',
                  './folder/s2_S011_L001_R1.fastq.gz',
                  './folder/s1_S009_L001_R1.fastq.gz']

        rev_fp = ['./folder/s4_S013_L001_R2.fastq.gz',
                  './folder/s2_S011_L001_R2.fastq.gz',
                  './folder/s1_S009_L001_R2.fastq.gz']

        with self.assertRaises(ValueError):
            make_read_pairs_per_sample(fwd_fp, rev_fp, fp)

    def test_make_read_pairs_per_sample_match_fwd_rev_notmatch(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)

        fwd_fp = ['./folder/s3_S013_L001_R1.fastq.gz',
                  './folder/s2_S011_L001_R1.fastq.gz',
                  './folder/s1_S009_L001_R1.fastq.gz']

        rev_fp = ['./folder/s3_S013_L001_R2.fastq.gz',
                  './folder/s2_S011_L001_R2.fastq.gz']

        with self.assertRaises(ValueError):
            make_read_pairs_per_sample(fwd_fp, rev_fp, fp)

    def test_make_read_pairs_per_sample_match_fwd_no_rp(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)

        fwd_fp = ['./folder/s3_S013_L001_R1.fastq.gz',
                  './folder/s2_S011_L001_R1.fastq.gz',
                  './folder/s4_S009_L001_R1.fastq.gz']

        rev_fp = []

        with self.assertRaises(ValueError):
            make_read_pairs_per_sample(fwd_fp, rev_fp, fp)

    def test_make_read_pairs_per_sample_match_fwd_2match(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)

        fwd_fp = ['./folder/s3_S013_L001_R1.fastq.gz',
                  './folder/s2_S011_L001_R1.fastq.gz',
                  './folder/s2_S009_L001_R1.fastq.gz']

        rev_fp = []

        with self.assertRaises(ValueError):
            make_read_pairs_per_sample(fwd_fp, rev_fp, fp)

    def test_generate_trim_analysis_commands_forward_reverse(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)

        exp_cmd = [
            'atropos trim -A GATCGGAAGAGCGTCGTGTAGGGAAAGGAGTGT '
            '--adapter GATCGGAAGAGCACACGTCTGAACTCCAGTCAC --max-n 80 '
            '--minimum-length 80 --pair-filter any --quality-cutoff 15 '
            '--threads 5 --trim-n -o output/s1.R1.fastq.gz '
            '-p output/s1.R2.fastq.gz -pe1 fastq/s1.fastq.gz '
            '-pe2 fastq/s1.R2.fastq.gz',
            'atropos trim -A GATCGGAAGAGCGTCGTGTAGGGAAAGGAGTGT '
            '--adapter GATCGGAAGAGCACACGTCTGAACTCCAGTCAC --max-n 80 '
            '--minimum-length 80 --pair-filter any --quality-cutoff 15 '
            '--threads 5 --trim-n -o output/s2.R1.fastq.gz '
            '-p output/s2.R2.fastq.gz -pe1 fastq/s2.fastq.gz '
            '-pe2 fastq/s2.R2.fastq.gz',
            'atropos trim -A GATCGGAAGAGCGTCGTGTAGGGAAAGGAGTGT '
            '--adapter GATCGGAAGAGCACACGTCTGAACTCCAGTCAC --max-n 80 '
            '--minimum-length 80 --pair-filter any --quality-cutoff 15 '
            '--threads 5 --trim-n -o output/s3.R1.fastq.gz '
            '-p output/s3.R2.fastq.gz -pe1 fastq/s3.fastq.gz '
            '-pe2 fastq/s3.R2.fastq.gz'
            ]

        exp_sample = [
            ('s1', 'SKB8.640193', 'fastq/s1.fastq.gz', 'fastq/s1.R2.fastq.gz'),
            ('s2', 'SKD8.640184', 'fastq/s2.fastq.gz', 'fastq/s2.R2.fastq.gz'),
            ('s3', 'SKB7.640196', 'fastq/s3.fastq.gz', 'fastq/s3.R2.fastq.gz')]

        obs_cmd, obs_sample = generate_trim_commands(
            ['fastq/s1.fastq.gz', 'fastq/s2.fastq.gz', 'fastq/s3.fastq.gz'],
            ['fastq/s1.R2.fastq.gz', 'fastq/s2.R2.fastq.gz',
             'fastq/s3.R2.fastq.gz'],
            fp, 'output', self.params)

        self.assertEqual(obs_cmd, exp_cmd)
        self.assertEqual(obs_sample, exp_sample)

    def test_trim(self):
        # generating filepaths
        in_dir = mkdtemp()
        self._clean_up_files.append(in_dir)

        fp1_1 = join(in_dir, 'kd_test_1_R1.fastq.gz')
        fp1_2 = join(in_dir, 'kd_test_1_R2.fastq.gz')
        fp2_1 = join(in_dir, 'kd_test_2_R1.fastq.gz')
        fp2_2 = join(in_dir, 'kd_test_2_R2.fastq.gz')
        copyfile('support_files/kd_test_1_R1.fastq.gz', fp1_1)
        copyfile('support_files/kd_test_1_R2.fastq.gz', fp1_2)
        copyfile('support_files/kd_test_2_R1.fastq.gz', fp2_1)
        copyfile('support_files/kd_test_2_R2.fastq.gz', fp2_2)

        # inserting new prep template
        prep_info_dict = {
            'SKB7.640196': {'run_prefix': 'kd_test_1'},
            'SKB8.640193': {'run_prefix': 'kd_test_2'}
        }
        data = {'prep_info': dumps(prep_info_dict),
                # magic #1 = testing study
                'study': 1,
                'data_type': 'Metagenomic'}
        pid = self.qclient.post('/apitest/prep_template/', data=data)['prep']

        # inserting artifacts
        data = {
            'filepaths': dumps([
                (fp1_1, 'raw_forward_seqs'),
                (fp1_2, 'raw_reverse_seqs'),
                (fp2_1, 'raw_forward_seqs'),
                (fp2_2, 'raw_reverse_seqs')]),
            'type': "per_sample_FASTQ",
            'name': "Test QC_Trim artifact",
            'prep': pid}
        aid = self.qclient.post('/apitest/artifact/', data=data)['artifact']

        self.params['input'] = aid
        data = {'user': 'demo@microbio.me',
                'command': dumps(['qp-woltka', '2020.11', 'Atropos v1.1.24']),
                'status': 'running',
                'parameters': dumps(self.params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        success, ainfo, msg = trim(self.qclient, jid, self.params, out_dir)

        self.assertEqual("", msg)
        self.assertTrue(success)

        # we are expecting 1 artifact in total
        self.assertEqual(1, len(ainfo))

        obs_fps = []
        for a in ainfo:
            self.assertEqual("per_sample_FASTQ", a.artifact_type)
            obs_fps.append(a.files)
        od = partial(join, out_dir)

        # ftype = 'per_sample_FASTQ'
        exp_fps = [
            [(od('kd_test_1.R1.fastq.gz'), 'raw_forward_seqs'),
             (od('kd_test_1.R2.fastq.gz'), 'raw_reverse_seqs'),
             (od('kd_test_2.R1.fastq.gz'), 'raw_forward_seqs'),
             (od('kd_test_2.R2.fastq.gz'), 'raw_reverse_seqs')]]
        self.assertEqual(exp_fps, obs_fps)

    def test_trim_just_fwd(self):
        # generating filepaths
        in_dir = mkdtemp()
        self._clean_up_files.append(in_dir)

        fp1 = join(in_dir, 'kd_test_1_R1.fastq.gz')
        fp2 = join(in_dir, 'kd_test_2_R1.fastq.gz')
        copyfile('support_files/kd_test_1_R1.fastq.gz', fp1)
        copyfile('support_files/kd_test_2_R1.fastq.gz', fp2)

        # inserting new prep template
        prep_info_dict = {
            'SKB7.640196': {'run_prefix': 'kd_test_1'},
            'SKB8.640193': {'run_prefix': 'kd_test_2'}
        }
        data = {'prep_info': dumps(prep_info_dict),
                # magic #1 = testing study
                'study': 1,
                'data_type': 'Metagenomic'}
        pid = self.qclient.post('/apitest/prep_template/', data=data)['prep']

        # inserting artifacts
        data = {
            'filepaths': dumps([
                (fp1, 'raw_forward_seqs'),
                (fp2, 'raw_forward_seqs')]),
            'type': "per_sample_FASTQ",
            'name': "Test QC_Trim artifact",
            'prep': pid}
        aid = self.qclient.post('/apitest/artifact/', data=data)['artifact']

        self.params['input'] = aid
        data = {'user': 'demo@microbio.me',
                'command': dumps(['qp-woltka', '2020.11', 'Atropos v1.1.24']),
                'status': 'running',
                'parameters': dumps(self.params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        success, ainfo, msg = trim(self.qclient, jid, self.params, out_dir)

        self.assertEqual("", msg)
        self.assertTrue(success)

        # we are expecting 1 artifact in total
        self.assertEqual(1, len(ainfo))

        obs_fps = []
        for a in ainfo:
            self.assertEqual("per_sample_FASTQ", a.artifact_type)
            obs_fps.append(a.files)
        od = partial(join, out_dir)

        # ftype = 'per_sample_FASTQ'
        exp_fps = [
            [(od('kd_test_1.R1.fastq.gz'), 'raw_forward_seqs'),
             (od('kd_test_2.R1.fastq.gz'), 'raw_forward_seqs')]]
        self.assertEqual(exp_fps, obs_fps)

    def test_per_sample_ainfo_error(self):
        in_dir = mkdtemp()
        self._clean_up_files.append(in_dir)
        makedirs(join(in_dir, 'sampleA'))
        makedirs(join(in_dir, 'sampleB'))

        # Paired-end
        with self.assertRaises(ValueError):
            _per_sample_ainfo(in_dir, (('sampleA', None, None, None),
                                       ('sampleB', None, None, None)), [],
                              'Atropos', 'QC_Trim Files', True)


MAPPING_FILE = (
    "#SampleID\tplatform\tbarcode\texperiment_design_description\t"
    "library_construction_protocol\tcenter_name\tprimer\trun_prefix\t"
    "instrument_model\tDescription\n"
    "SKB7.640196\tILLUMINA\tA\tA\tA\tANL\tA\ts3\tIllumina MiSeq\tdesc1\n"
    "SKB8.640193\tILLUMINA\tA\tA\tA\tANL\tA\ts1\tIllumina MiSeq\tdesc2\n"
    "SKD8.640184\tILLUMINA\tA\tA\tA\tANL\tA\ts2\tIllumina MiSeq\tdesc3\n"
)

MAPPING_FILE_2 = (
    "#SampleID\tplatform\tbarcode\texperiment_design_description\t"
    "library_construction_protocol\tcenter_name\tprimer\t"
    "run_prefix\tinstrument_model\tDescription\n"
    "SKB7.640196\tILLUMINA\tA\tA\tA\tANL\tA\ts3\tIllumina MiSeq\tdesc1\n"
    "SKB8.640193\tILLUMINA\tA\tA\tA\tANL\tA\ts1\tIllumina MiSeq\tdesc2\n"
    "SKD8.640184\tILLUMINA\tA\tA\tA\tANL\tA\ts1\tIllumina MiSeq\tdesc3\n"
)


if __name__ == '__main__':
    main()
