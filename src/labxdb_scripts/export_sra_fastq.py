#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# Copyright © 2016 Charles E. Vejnar
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
#

import argparse
import csv
import json
import os
import subprocess
import sys
import time

import labxdb
import labxdb.fastq

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
    parser.add_argument('--fastq_exts', dest='fastq_exts', action='store', default='.fastq,.fastq.gz,.fastq.zst', help='FASTQ file extensions (comma separated).')
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
        if a in ['fastq_exts']:
            if v is None:
                config[a] = []
            else:
                config[a] = [r.strip() for r in v.split(',')]

    # Open data input info
    data = json.load(open(config['path_data']))

    # Prepare jobs
    jobs = []
    for record in data:
        if record['type'] == 'single':
            # List FASTQ files
            fastqs = labxdb.fastq.find_fastqs(os.path.join(config['path_seq_run'], record['run_ref']), config['fastq_exts'], [labxdb.fastq.parse_illumina_fastq_filename, labxdb.fastq.parse_fastq_filename, labxdb.fastq.parse_no_end_fastq_filename])
            # Create job(s)
            for name, path_list in fastqs.items():
                for p in path_list:
                    path_in = os.path.join(p['path'], p['fname'])
                    if p['end'] is None:
                        fname_out = f"{record['run_ref']}.fastq"
                    else:
                        fname_out = f"{record['run_ref']}_{p['end']}.fastq"
                    # Add zip job
                    if config['overwrite'] or (config['overwrite'] is False and os.path.exists(fname_out+'.gz') is False):
                        jobs.append([path_in, fname_out])
                    # Add file to data
                    if 'files' in record:
                        record['files'].append(fname_out)
                    else:
                        record['files'] = [fname_out]
        elif record['type'] == 'multiplex':
            # List FASTQ files
            fastqs = labxdb.fastq.find_fastqs(os.path.join(config['path_seq_raw'], record['flowcell']), config['fastq_exts'], [labxdb.fastq.parse_illumina_fastq_filename, labxdb.fastq.parse_fastq_filename, labxdb.fastq.parse_no_end_fastq_filename])
            # Create job(s)
            for name, path_list in fastqs.items():
                if name == record['tube_label']:
                    for p in path_list:
                        path_in = os.path.join(p['path'], p['fname'])
                        if p['end'] is None:
                            fname_out = f"{record['flowcell']}_{record['tube_label']}.fastq"
                        else:
                            fname_out = f"{record['flowcell']}_{record['tube_label']}_{p['end']}.fastq"
                        # Add zip job
                        if config['overwrite'] or (config['overwrite'] is False and os.path.exists(fname_out+'.gz') is False):
                            jobs.append([path_in, fname_out])
                        # Add file to data
                        if 'files' in record:
                            record['files'].append(fname_out)
                        else:
                            record['files'] = [fname_out]

    # Replace filename fields in data TSV
    # Open CSV
    csv.register_dialect('excel-tab-unix', delimiter='\t', lineterminator='\n', quoting=csv.QUOTE_MINIMAL)
    rcsv = csv.reader(open('sra_data.tsv', 'rt'), dialect='excel-tab-unix')
    csv_headers = next(rcsv)
    # Max number of files per run
    max_file = max([len(record['files']) for record in data])
    # Filename headers
    csv_headers_filename = [i for i in range(len(csv_headers)) if csv_headers[i].startswith('filename')]
    assert max_file <= len(csv_headers_filename), f'Insufficient ({max_file}) filename column(s) found ({len(csv_headers_filename)})'
    # Add files
    records = []
    for irecord, record in enumerate(rcsv):
        for i in range(len(data[irecord]['files'])):
            assert data[irecord]['run_ref'] == record[csv_headers.index('library_ID')]
            record[csv_headers_filename[i]] = data[irecord]['files'][i] + '.gz'
        records.append(record)
    # Save old TSV
    os.replace('sra_data.tsv', 'sra_data.tsv.backup')
    # Write TSV
    with open('sra_data.tsv', 'w', newline='') as f:
        writer = csv.writer(f, dialect='excel-tab-unix')
        writer.writerow(csv_headers)
        writer.writerows(records)

    # Run jobs
    pfu.parallel.run(rezip, jobs, config['num_processor'])

if __name__ == '__main__':
    sys.exit(main())
