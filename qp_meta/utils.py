# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# This file contains functions used by multiple commands
# -----------------------------------------------------------------------------
from qiita_client.util import system_call, get_sample_names_by_run_prefix
from itertools import zip_longest
from os.path import basename, join, exists
from functools import partial
from qiita_client import ArtifactInfo


def make_read_pairs_per_sample(forward_seqs, reverse_seqs, map_file):
    """Recovers read pairing information

    Parameters
    ----------
    forward_seqs : list of strs
        The list of forward seqs filepaths
    reverse_seqs : list of str
        The list of reverse seqs filepaths
    map_file : str
        The path to the mapping file

    Returns
    -------
    samples: list of tup
        list of 4-tuples with run prefix, sample name, fwd read fp, rev read fp

    Raises
    ------
    ValueError
        If the rev is not an empty list and the same length than fwd seqs
        The prefixes of the run_prefix don't match the file names

    Notes
    -----
    At this stage it is required that if reverse sequences are present that all
    samples have both a forward and a reverse sequence. However, the read
    trimming/filtering step can sometimes eliminate all reverse reads,
    especially in low
    coverage samples with poor overall reverse read quality.
    """

    # sort forward seqs
    forward_seqs.sort()

    # check that rev seqs are same len
    if reverse_seqs:
        if len(forward_seqs) != len(reverse_seqs):
            raise ValueError('Your reverse and forward files are of different '
                             'length. Forward: %s. Reverse: %s.' %
                             (', '.join(forward_seqs),
                              ', '.join(reverse_seqs)))
        reverse_seqs.sort()

    # get run prefixes
    # These are prefixes that should match uniquely to forward reads
    # sn_by_rp is dict of samples keyed by run prefixes
    sn_by_rp = get_sample_names_by_run_prefix(map_file)

    # make pairings
    samples = []
    used_prefixes = set()
    for i, (fwd_fp, rev_fp) in enumerate(zip_longest(forward_seqs,
                                                     reverse_seqs)):
        # fwd_fp is the fwd read filepath
        fwd_fn = basename(fwd_fp)

        # iterate over run prefixes and make sure only one matches
        run_prefix = None
        for rp in sn_by_rp:
            if fwd_fn.startswith(rp) and run_prefix is None:
                run_prefix = rp
            elif fwd_fn.startswith(rp) and run_prefix is not None:
                raise ValueError('Multiple run prefixes match this fwd read: '
                                 '%s' % fwd_fn)

        # make sure that we got one matching run prefix:
        if run_prefix is None:
            raise ValueError('No run prefix matching this fwd read: %s'
                             % fwd_fn)

        if run_prefix in used_prefixes:
            raise ValueError('This run prefix matches multiple fwd reads: '
                             '%s' % run_prefix)

        if rev_fp is None:
            samples.append((run_prefix, sn_by_rp[run_prefix], fwd_fp, None))
        else:
            rev_fn = basename(rev_fp)
            # if we have reverse reads, make sure the matching pair also
            # matches the run prefix:
            if not rev_fn.startswith(run_prefix):
                raise ValueError('Reverse read does not match this run prefix.'
                                 '\nRun prefix: %s\nForward read: %s\n'
                                 'Reverse read: %s\n' %
                                 (run_prefix, fwd_fn, rev_fn))

            samples.append((run_prefix, sn_by_rp[run_prefix], fwd_fp,
                            rev_fp))

        used_prefixes.add(run_prefix)

    return(samples)


def _format_params(parameters, func_params):
    params = []
    # Loop through all of the commands alphabetically
    for param in sorted(func_params):
        # Find the value using long parameter names
        parameter = func_params[param]
        value = parameters[parameter]
        dash = '--'
        # Check for single letter commands
        if len(param) == 1:
            dash = '-'
        if value is True:
            params.append('%s%s' % (dash, param))
        elif value is False:
            continue
        elif value and value != 'default':
            params.append('%s%s %s' % (dash, param, value))

    param_string = ' '.join(params)

    return(param_string)


def _run_commands(qclient, job_id, commands, msg, cmd_name):
    for i, cmd in enumerate(commands):
        qclient.update_job_step(job_id, msg % (i+1))
        std_out, std_err, return_value = system_call(cmd)
        if return_value != 0:
            error_msg = ("Error running %s:\nStd out: %s\nStd err: %s"
                         "\n\nCommand run was:\n%s" % (cmd_name, std_out,
                                                       std_err, cmd))
            return False, error_msg

    return True, ""


def _per_sample_ainfo(
        out_dir, samples, suffixes, prg_name,
        files_type_name, fwd_and_rev=False):
    files = []
    missing_files = []
    smd = partial(join, out_dir)
    for rp, _, _, _ in samples:
        for suff in suffixes:
            fname = smd(suff % rp)
            if exists(fname):
                if fname.endswith('R1.fastq.gz'):
                    ftype = 'raw_forward_seqs'
                elif fname.endswith('R2.fastq.gz'):
                    ftype = 'raw_reverse_seqs'
                else:
                    # this should never happen and it's not really possible
                    # to reproduce so no tests!
                    raise ValueError('File %s has an unexpected name' % fname)
                files.append((fname, ftype))
            else:
                missing_files.append(fname)

    if not files:
        # Command did not create any files, which means that no sequence
        # was kept after quality control and filtering for host data
        raise ValueError("No sequences left after %s" % prg_name)

    return [ArtifactInfo(files_type_name, 'per_sample_FASTQ', files)]
