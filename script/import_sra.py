#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# Copyright (C) 2015-2020 Charles E. Vejnar
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
#

import argparse
import glob
import json
import os
import re
import shutil
import subprocess
import sys

import labxdb
import labxdb.ncbi

import pyfnutils as pfu
import pyfnutils.parallel

def print_summary(samples):
    for sample in samples:
        print('{:<20}{:<15}'.format(*[str(sample[s]) for s in ['ref', 'label']]))
        for run in sample['runs']:
            print('   \_'+run['ref'])
    print()

def add_samples(samples, srp, srp_title, dbl):
    # Add project
    dbl.post('project/new', json=[{'project_ref':srp, 'label_short':srp, 'label_long':srp_title, 'sra_ref':srp}])

    # Add samples
    query = []
    for sample in samples:
        # Add run(s)
        for irun, run in enumerate(sample['runs']):
            # Add run
            runs = dbl.post('run/new', json=[{'run_ref':run['ref'], 'run_order':irun+1, 'failed':False, 'platform':run['platform'], 'paired':run['paired'], 'sra_ref':run['ref']}])
            # Save new run ID
            run['id'] = runs[0]['run_id']
            # Add run to query
            query.append([['new', run['ref']], ['new', sample['label']], ['new', sample['label']], ['append', srp]])
    result = dbl.post('assign', json={'prefix':'SI', 'new_run':False, 'query':query})

    # Update replicate sra_ref
    for ireplicate, replicate in enumerate(result['ids'][1]):
        dbl.post('replicate/edit/'+str(replicate), json=[{'sra_ref':samples[ireplicate]['ref']}])

def dump_sra(samples, config, dbl):
    # Create output dir
    path_out = os.path.join(config['path_seq_prepared'], config['project'])
    if not os.path.exists(path_out):
        os.mkdir(path_out)
    # Dump and split SRA file(s)
    jobs = [[run['ref'], run.get('id'), run.get('paired'), path_out, config, dbl] for sample in samples for run in sample['runs']]
    pfu.parallel.run(dump_sra_run, jobs, num_processor=config['num_processor'])

def dump_sra_run(srr, run_id, paired, path_out, config, dbl):
    # Prepare command
    cmd = ['fasterq-dump', '--progress']
    if paired is True:
        cmd.append('--split-3')
    cmd.extend(['--outdir', path_out, '--temp', path_out, '--threads', str(config['num_processor'])])
    cmd.append(srr)
    # Run command
    print('Download and dump %s'%srr)
    subprocess.run(cmd, check=True)
    # Delete orphan reads or Rename single end
    if os.path.exists(os.path.join(path_out, srr+'_1.fastq')):
        path_orphan = os.path.join(path_out, srr+'.fastq')
        if os.path.exists(path_orphan):
            print('Remove', path_orphan)
            os.remove(path_orphan)
    else:
        os.rename(os.path.join(path_out, srr+'.fastq'), os.path.join(path_out, srr+'_1.fastq'))
    # Post-dump
    for ifq, fq in enumerate(sorted(glob.glob(os.path.join(path_out, srr+'_*.fastq')))):
        fq = os.path.basename(fq)
        prefix, readn = re.search(r'(.*)_(\d+).fastq', fq).groups()
        fn = '%s_R%s.fastq'%(prefix, readn)
        # Strip
        p = subprocess.run(['fastq_strip', fq, fn], check=True, capture_output=True, text=True, cwd=path_out)
        nread, max_length = p.stdout.strip().split()
        # Update max length in database
        if run_id is None:
            run = dbl.get('run/get-ref/'+srr)
            if len(run[0]) > 0 and len(run[0][0]) > 0:
                run_id = run[0][0]['run_id']
        if run_id is not None and nread is not None and max_length is not None:
            dbl.post('run/edit/'+str(run_id), json=[{'spots':int(nread), 'max_read_length':int(max_length)}])
        # Remove
        os.remove(os.path.join(path_out, fq))
        # Zip
        if 'zip_cmd' in config:
            subprocess.run(config['zip_cmd'] + [fn], check=True, cwd=path_out)
    # Read-only permission
    for f in glob.glob(os.path.join(path_out, srr+'*.fastq*')):
        os.chmod(os.path.join(path_out, os.path.basename(f)), 0o0444)

def create_links(samples, config):
    path_out = os.path.join(config['path_seq_prepared'], config['project'])
    # Link to FASTQ file(s)
    for sample in samples:
        for run in sample['runs']:
            srr = run['ref']
            print('Link %s'%srr)
            # Create run folder
            path_run = os.path.join(config['path_seq_run'], srr)
            os.makedirs(path_run, exist_ok=True)
            # Create links
            for fq in glob.glob(os.path.join(path_out, '%s*'%srr)):
                fq_rel_path = os.path.relpath(fq, path_run)
                fq_new_path = os.path.join(path_run, os.path.basename(fq))
                if not os.path.exists(fq_new_path) or config['overwrite_links']:
                    os.symlink(fq_rel_path, fq_new_path)

def check_exe(names):
    for name in names:
        if shutil.which(name) == None:
            raise FileNotFoundError('%s missing'%name)

def main(argv=None):
    if argv is None:
        argv = sys.argv
    parser = argparse.ArgumentParser(description='Import SRA data.')
    parser.add_argument('-j', '--project', dest='project', action='store', required=True, help='SRPxxx.')
    parser.add_argument('-i', '--db_import', dest='db_import', action='store_true', help='Import into database.')
    parser.add_argument('-d', '--dump_sra', dest='dump_sra', action='store_true', help='Download and dump SRA file(s) to FASTQ.')
    parser.add_argument('-z', '--zip_cmd', dest='zip_cmd', action='store', default='zstd,--rm,-19,-T4', help='Zip command (comma separated).')
    parser.add_argument('-l', '--create_links', dest='create_links', action='store_true', help='Create link(s) to FASTQ.')
    parser.add_argument('-a', '--path_seq_prepared', dest='path_seq_prepared', action='store', help='Path to prepared.')
    parser.add_argument('-n', '--path_seq_run', dest='path_seq_run', action='store', help='Path to run.')
    parser.add_argument('-r', '--runs', dest='runs', action='store', help='Selected runs excluding all others in the project (comma separated).')
    parser.add_argument('-s', '--save_sra_xml', dest='save_sra_xml', action='store_true', help='Save project XML from SRA.')
    parser.add_argument('-p', '--processor', dest='num_processor', action='store', type=int, default=1, help='Number of processor')
    parser.add_argument('--overwrite_links', dest='overwrite_links', action='store_true', default=False, help='Overwrite existing links (default:False).')
    parser.add_argument('--path_config', dest='path_config', action='store', help='Path to config')
    parser.add_argument('--http_url', '--labxdb_http_url', dest='labxdb_http_url', action='store', help='Database HTTP URL')
    parser.add_argument('--http_login', '--labxdb_http_login', dest='labxdb_http_login', action='store', help='Database HTTP login')
    parser.add_argument('--http_password', '--labxdb_http_password', dest='labxdb_http_password', action='store', help='Database HTTP password')
    parser.add_argument('--http_path', '--labxdb_http_path', dest='labxdb_http_path', action='store', help='Database HTTP path')
    parser.add_argument('--http_db', '--labxdb_http_db', dest='labxdb_http_db', action='store', help='Database HTTP DB')
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
        if a in ['zip_cmd', 'runs']:
            if v is None:
                config[a] = []
            else:
                config[a] = [r.strip() for r in v.split(',')]

    # Checks
    if config['dump_sra']:
        # Check folder
        if not os.path.exists(config['path_seq_prepared']):
            print('ERROR: %s not found'%config['path_seq_prepared'])
            return 1
        # Check executables
        exes = ['fasterq-dump', 'fastq_strip']
        if 'zip_cmd' in config:
            exes.append(config['zip_cmd'][0])
        check_exe(exes)

    # Init. DBLink
    if 'labxdb_http_path' not in config and 'labxdb_http_db' not in config:
        if 'labxdb_http_path_seq' in config:
            config['labxdb_http_path'] = config['labxdb_http_path_seq']
        else:
            config['labxdb_http_db'] = 'seq'
    dbl = labxdb.DBLink(config.get('labxdb_http_url'), config.get('labxdb_http_login'), config.get('labxdb_http_password'), config.get('labxdb_http_path'), config.get('labxdb_http_db'))

    # Get info
    srp_title, samples = labxdb.ncbi.get_samples_infos(config['project'], import_runs=config['runs'], save_sra_xml=config['save_sra_xml'], verbose=True)
    print('\n'+srp_title)
    print_summary(samples)

    # Import
    if config['db_import']:
        print('Import to DB')
        add_samples(samples, config['project'], srp_title, dbl)
    if config['dump_sra']:
        print('Dump from SRA')
        dump_sra(samples, config, dbl)
    if config['create_links']:
        print('Create links')
        create_links(samples, config)

if __name__ == '__main__':
    sys.exit(main())
