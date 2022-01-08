#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# Copyright (C) 2016-2022 Charles E. Vejnar
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
#

import argparse
import glob
import json
import os
import subprocess
import sys
import time

import pyfnutils as pfu
import pyfnutils.parallel

def rezip(path_in, fname_out):
    print(fname_out)
    subprocess.run(['zstd', '-dkf', path_in, '-o', fname_out], check=True)
    time.sleep(10)
    subprocess.run(['gzip', '-f9', fname_out], check=True)

def main(argv=None):
    if argv is None:
        argv = sys.argv
    parser = argparse.ArgumentParser(description='Export data to SRA (prepare FASTQ files).')
    parser.add_argument('-r', '--path_seq_raw', dest='path_seq_raw', action='store', help='Path to raw root')
    parser.add_argument('-s', '--path_seq_run', dest='path_seq_run', action='store', help='Path to run root')
    parser.add_argument('-d', '--path_data', dest='path_data', action='store', default='data.json', help='Run data')
    parser.add_argument('-c', '--path_config', dest='path_config', action='store', help='Path to config')
    parser.add_argument('-f', '--overwrite', dest='overwrite', action='store_true', help='Overwrite already exported files')
    parser.add_argument('-p', '--processor', dest='num_processor', action='store', type=int, default=2, help='Number of processor')
    args = parser.parse_args(argv[1:])

    # Get config (JSON single file or all files in path_config)
    config = {}
    paths = []
    if args.path_config is None:
        if 'HTS_CONFIG_PATH' in os.environ:
            paths.append(os.environ['HTS_CONFIG_PATH'])
        elif 'XDG_CONFIG_HOME' in os.environ:
            paths.append(os.path.join(os.environ['XDG_CONFIG_HOME'], 'hts'))
    else:
        paths.append(args.path_config)
    for path in paths:
        if os.path.isdir(path):
            for f in sorted(os.listdir(path)):
                if f.endswith('.json'):
                    config = {**config, **json.load(open(os.path.join(path, f)))}
        elif os.path.isfile(path):
            config = {**config, **json.load(open(path))}

    # Input local config from args
    vargs = vars(args)
    for a, v in vargs.items():
        if v is not None and (a not in config or v != parser.get_default(a)):
            config[a] = v

    # Open data input info
    data = json.load(open(config['path_data']))

    # Prepare jobs
    jobs = []
    for record in data:
        if record['type'] == 'single':
            gin = glob.glob(os.path.join(config['path_seq_run'], record['run_ref'], f"*{record['pattern']}.zst"))
            if len(gin) == 1:
                path_in = gin[0]
            else:
                raise ValueError('More than one file found')
            fname_out = f"{record['run_ref']}{record['pattern']}"
            if config['overwrite'] or (config['overwrite'] is False and os.path.exists(fname_out+'.gz') is False):
                jobs.append([path_in, fname_out])
        elif record['type'] == 'multiplex':
            gin = glob.glob(os.path.join(config['path_seq_raw'], record['flowcell'], f"{record['tube_label']}*{record['pattern']}.zst"))
            if len(gin) == 1:
                path_in = gin[0]
            else:
                raise ValueError('More than one file found')
            fname_out = f"{record['flowcell']}_{record['tube_label']}{record['pattern']}"
            if config['overwrite'] or (config['overwrite'] is False and os.path.exists(fname_out+'.gz') is False):
                jobs.append([path_in, fname_out])

    # Run jobs
    pfu.parallel.run(rezip, jobs, config['num_processor'])

if __name__ == '__main__':
    sys.exit(main())
