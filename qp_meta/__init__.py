# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from qiita_client import QiitaPlugin

from .sortmerna import sortmerna_cmd
from .utils import plugin_details

# Initialize the plugin
plugin = QiitaPlugin(**plugin_details)

plugin.register_command(sortmerna_cmd)
