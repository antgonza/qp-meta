# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from qiita_client import QiitaPlugin

from .sortmerna import sortmerna_cmd


# Initialize the plugin
plugin = QiitaPlugin(
    'qp-meta', '2021.01', 'meta analysis tools for shotgun data')

plugin.register_command(sortmerna_cmd)
