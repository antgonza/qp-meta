# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from qiita_client import QiitaPlugin

from .trim import trim_cmd
from .filter import filter_cmd
from .woltka import woltka_cmd
from .sortmerna import sortmerna_cmd


# Initialize the plugin
plugin = QiitaPlugin(
    'qp-woltka', '2020.11', 'woltka analysis tools for shotgun data')

plugin.register_command(trim_cmd)
plugin.register_command(filter_cmd)
plugin.register_command(woltka_cmd)
plugin.register_command(sortmerna_cmd)
