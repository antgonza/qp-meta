# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import main
from os import close, remove, makedirs
from os.path import exists, isdir, join
from shutil import rmtree, copyfile
from tempfile import mkstemp, mkdtemp
from json import dumps
from functools import partial
import os

from qiita_client.testing import PluginTestCase

from qp_meta import plugin
from qp_meta.sortmerna.sortmerna import (
    generate_sortmerna_commands, sortmerna)
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

    def test_generate_sortmerna_analysis_commands_forward_reverse(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)

        db_path = os.environ["QC_SORTMERNA_DB_DP"]

        rna_ref_db = (
                '{0}silva-arc-23s-id98.fasta,'
                '{0}silva-arc-23s-id98.idx:'
                '{0}silva-bac-16s-id90.fasta,'
                '{0}silva-bac-16s-id90.idx:'
                '{0}silva-bac-23s-id98.fasta,'
                '{0}silva-bac-23s-id98.idx:'
                '{0}silva-arc-16s-id95.fasta,'
                '{0}silva-arc-16s-id95.idx:'
                '{0}silva-euk-18s-id95.fasta,'
                '{0}silva-euk-18s-id95.idx:'
                '{0}silva-euk-28s-id98.fasta,'
                '{0}silva-euk-28s-id98.idx:'
                '{0}rfam-5s-database-id98.fasta,'
                '{0}rfam-5s-database-id98.idx:'
                '{0}rfam-5.8s-database-id98.fasta,'
                '{0}rfam-5.8s-database-id98.idx'
        ).format(db_path)
        # index files take up the most space

        exp_cmd = [
            ('unpigz -p 5 -c fastq/s1.fastq.gz > output/s1.fastq && '

             'sortmerna --ref %s --reads output/s1.fastq '
             '--aligned output/s1.ribosomal.R1 '
             '--other output/s1.nonribosomal.R1 '
             '--fastx -a 5 --blast 1 -m 3988 --num_alignments 1 && '

             'pigz -p 5 -c output/s1.ribosomal.R1.fastq > '
             'output/s1.ribosomal.R1.fastq.gz && '

             'pigz -p 5 -c output/s1.nonribosomal.R1.fastq > '
             'output/s1.nonribosomal.R1.fastq.gz;') % rna_ref_db,
            ('unpigz -p 5 -c fastq/s1.R2.fastq.gz > output/s1.R2.fastq && '

             'sortmerna --ref %s --reads output/s1.R2.fastq '
             '--aligned output/s1.ribosomal.R2 '
             '--other output/s1.nonribosomal.R2 '
             '--fastx -a 5 --blast 1 -m 3988 --num_alignments 1 && '

             'pigz -p 5 -c output/s1.ribosomal.R2.fastq > '
             'output/s1.ribosomal.R2.fastq.gz && '

             'pigz -p 5 -c output/s1.nonribosomal.R2.fastq > '
             'output/s1.nonribosomal.R2.fastq.gz;') % rna_ref_db
        ]

        exp_sample = [
            ('s1', 'SKB8.640193', 'fastq/s1.fastq.gz', 'fastq/s1.R2.fastq.gz')
            ]
        obs_cmd, obs_sample = generate_sortmerna_commands(
            ['fastq/s1.fastq.gz'], ['fastq/s1.R2.fastq.gz'],
            fp, 'output', self.params)

        self.assertEqual(obs_cmd, exp_cmd)
        self.assertEqual(obs_sample, exp_sample)

    def test_sortmerna(self):
        # generating filepaths
        in_dir = mkdtemp()
        self._clean_up_files.append(in_dir)

        fp1_1 = join(in_dir, 'kd_test_1_R1.fastq.gz')
        fp1_2 = join(in_dir, 'kd_test_1_R2.fastq.gz')
        copyfile('support_files/kd_test_1_R1.fastq.gz', fp1_1)
        copyfile('support_files/kd_test_1_R2.fastq.gz', fp1_2)
        # inserting new prep template
        prep_info_dict = {
            'SKB8.640193': {'run_prefix': 'kd_test_1'}
        }
        data = {'prep_info': dumps(prep_info_dict),
                'study': 1,
                'data_type': 'Metagenomic'}
        pid = self.qclient.post('/apitest/prep_template/', data=data)['prep']

        # inserting artifacts
        data = {
            'filepaths': dumps([
                (fp1_1, 'raw_forward_seqs'),
                (fp1_2, 'raw_reverse_seqs')]),
            'type': "per_sample_FASTQ",
            'name': "Test QC_Sortmerna artifact",
            'prep': pid}
        aid = self.qclient.post('/apitest/artifact/', data=data)['artifact']

        self.params['input'] = aid
        data = {'user': 'demo@microbio.me',
                'command': dumps(['qp-meta', '2021.01', 'Sortmerna v2.1b']),
                'status': 'running',
                'parameters': dumps(self.params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        success, ainfo, msg = sortmerna(self.qclient, jid, self.params,
                                        out_dir)

        self.assertEqual("", msg)
        self.assertTrue(success)

        # we are expecting 2 artifacts in total
        self.assertEqual(2, len(ainfo))
        obs_fps = []
        # checking if this is a dict of dicts
        for a in ainfo:
            self.assertEqual("per_sample_FASTQ", a.artifact_type)
            obs_fps.append(a.files)

        od = partial(join, out_dir)

        exp_fps = [
            [(od('kd_test_1.nonribosomal.R1.fastq.gz'), 'raw_forward_seqs'),
             (od('kd_test_1.nonribosomal.R2.fastq.gz'), 'raw_reverse_seqs')],
            [(od('kd_test_1.ribosomal.R1.fastq.gz'), 'raw_forward_seqs'),
             (od('kd_test_1.ribosomal.R2.fastq.gz'), 'raw_reverse_seqs')]]

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
                              'Sortmerna', 'QC_Sortmerna Files', True)


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
