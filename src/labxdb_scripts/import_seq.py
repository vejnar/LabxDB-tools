#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# Copyright (C) 2018-2022 Charles E. Vejnar
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
#

import argparse
import datetime
import glob
import json
import logging
import os
import shutil
import subprocess
import sys

import labxdb
import labxdb.fastq

import pyfnutils as pfu
import pyfnutils.log

class Error(Exception):
    def __init__(self, message):
        self.message = message

def download_read_files(url_remote, path_output, config, do_format=False, download_program='wget', no_readonly=False, dry_run=False, logger=None):
    # Parameters
    if logger is None:
        import logging as logger

    if download_program == 'wget':
        for ur in url_remote:
            cmds = ['wget', '--no-parent', '-m', ur]
            if 'wget_options' in config:
                cmds += config['wget_options']
            if dry_run:
                logger.info(f'DRY RUN: {cmds}')
            else:
                logger.info(f'Downloading with {cmds}')
                if not os.path.exists(path_output):
                    os.mkdir(path_output)
                try:
                    subprocess.run(cmds, check=True, cwd=path_output)
                except Exception as e:
                    if e.returncode == 8:
                        logger.warning('Some links were invalid')
                    else:
                        raise

    # Lock read-only output path
    if not dry_run and do_format is False and no_readonly is False:
        os.chmod(path_output, 0o0555)

def format_raw_read_files(bulk, path_seq_tmp, path_seq_raw, do_format='raw', delete_download=True, squashfs_download=False, flowcell_dir=True, fastq_exts=['.fastq'], no_readonly=False, dry_run=False, num_processor=1, logger=None):
    # Parameters
    if logger is None:
        import logging as logger

    # List FASTQ files
    if path_seq_tmp is not None:
        path_bulk = os.path.join(path_seq_tmp, bulk)
    else:
        path_bulk = os.path.join(path_seq_raw, bulk)
    fastqs = labxdb.fastq.find_fastqs(path_bulk, fastq_exts, [labxdb.fastq.parse_illumina_fastq_filename, labxdb.fastq.parse_fastq_filename])
    if labxdb.fastq.check_fastqs(fastqs) is False:
        raise Error('Name collision: Runs with the same name detected in different folder')

    # Reformat
    logger.info(f'Reformating following rule «{do_format}»')
    if do_format.startswith('run_'):
        for name, path_list in fastqs.items():
            path_list_per_end = {}
            for p in path_list:
                if p['end'] in path_list_per_end:
                    path_list_per_end[p['end']].append(os.path.join(p['path'], p['fname']))
                else:
                    path_list_per_end[p['end']] = [os.path.join(p['path'], p['fname'])]
            for end, paths in path_list_per_end.items():
                fname_out = os.path.join(path_bulk, f'{name}_{end}.fastq')
                # Concatenate
                cmd = ' '.join(['zcat'] + paths + ['>', fname_out])
                if dry_run:
                    logger.info(f'DRY RUN: Concatenating run with {cmd}')
                else:
                    logger.info(f'Concatenating run with {cmd}')
                    subprocess.run(cmd, check=True, shell=True)
                # Zip
                if do_format == 'run_zstd':
                    cmd = ['zstd', '--rm', '-T'+str(num_processor), '-19', fname_out]
                    if dry_run:
                        logger.info(f'DRY RUN: Zip run with {cmd}')
                    else:
                        logger.info(f'Zip run with {cmd}')
                        subprocess.run(cmd, check=True)
                        if no_readonly is False:
                            os.chmod(fname_out+'.zst', 0o0444)
                # Clean
                for p in paths:
                    if dry_run:
                        logger.info(f'DRY RUN: Removing {p}')
                    else:
                        logger.info(f'Removing {p}')
                        os.remove(p)

    # Squashfs
    path_download = os.path.join(path_bulk, 'download')
    path_squashfs = os.path.join(path_bulk, 'archive.sqfs')
    if squashfs_download:
        cmd = ['mksquashfs', '.', '../archive.sqfs', '-processors', str(num_processor), '-comp', 'xz', '-Xdict-size', '100%', '-b', '1M']
        if dry_run:
            logger.info(f'DRY RUN: Squashfs with {cmd}')
        else:
            logger.info(f'Squashfs with {cmd}')
            subprocess.run(cmd, check=True, cwd=path_download)
            # Set cleaning
            delete_download = True
            # Permission
            if no_readonly is False:
                os.chmod(path_squashfs, 0o0444)
    if delete_download:
        if dry_run:
            logger.info(f'DRY RUN: Cleaning {path_download}')
        else:
            # Stop logging into file
            for h in logger.handlers:
                if isinstance(h, logging.FileHandler):
                    h.close()
                    logger.removeHandler(h)
            # Remove download directory
            logger.info(f'Cleaning {path_download}')
            shutil.rmtree(path_download)

    # Update list of FASTQ files
    fastqs = labxdb.fastq.find_fastqs(path_bulk, fastq_exts, [labxdb.fastq.parse_illumina_fastq_filename, labxdb.fastq.parse_fastq_filename])

    # Use flowcell as bulk name
    if flowcell_dir:
        flowcells = set()
        for name, path_list in fastqs.items():
            # Extract flowcell name from first FASTQ
            r = labxdb.fastq.get_illumina_fastq_info(os.path.join(path_list[0]['path'], path_list[0]['fname']))
            flowcell = r['flowcell'].split(':')[2]
            flowcells.add(flowcell)
            # Create new flowcell directory
            path_new_bulk = os.path.join(path_seq_raw, flowcell)
            if dry_run:
                if not os.path.exists(path_new_bulk):
                    logger.info(f'DRY RUN: New directory for {flowcell}')
            else:
                if not os.path.exists(path_new_bulk):
                    logger.info(f'New directory for {flowcell}')
                    os.mkdir(path_new_bulk)
                elif not os.access(path_new_bulk, os.W_OK):
                    os.chmod(path_new_bulk, 0o0775)
            # Move archive.sqfs to new flowcell directory
            if os.path.exists(path_squashfs) and not os.path.exists(os.path.join(path_new_bulk, 'archive.sqfs')):
                shutil.move(path_squashfs, path_new_bulk)
            # Move file(s) to flowcell directory
            if dry_run:
                logger.info(f'DRY RUN: Move {name} from {bulk} to {path_new_bulk}')
            else:
                logger.info(f'Move {name} from {bulk} to {path_new_bulk}')
                for p in path_list:
                    shutil.move(os.path.join(path_bulk, p['fname']), path_new_bulk)
                    # Change path in fastqs too
                    p['path'] = path_new_bulk
                # Permission
                if no_readonly is False:
                    os.chmod(path_new_bulk, 0o0555)
        # Remove bulk directory
        if dry_run:
            logger.info(f'DRY RUN: Removing {path_bulk}')
        else:
            logger.info(f'Removing {path_bulk}')
            shutil.rmtree(path_bulk)
        if len(flowcells) > 1:
            logger.warning(f"Data from multiple flowcells: {','.join(flowcells)}")
    elif path_seq_tmp is not None:
        if dry_run:
            logger.info(f'DRY RUN: Move {path_bulk} to {path_seq_raw}')
        else:
            logger.info(f'Move {path_bulk} to {path_seq_raw}')
            shutil.move(path_bulk, path_seq_raw)
            # Permission
            if no_readonly is False:
                os.chmod(os.path.join(path_seq_raw, bulk), 0o0555)

    return fastqs

def import_staging(bulk, path_seq_raw, fastqs=None, ref_prefix='TMP_', fastq_exts=['.fastq'], dry_run=False, dbl=None, logger=None):
    # Parameters
    if logger is None:
        import logging as logger

    logger.info('Adding runs to table')
    if fastqs is None:
        fastqs = labxdb.fastq.find_fastqs(os.path.join(path_seq_raw, bulk), fastq_exts, [labxdb.fastq.parse_illumina_fastq_filename, labxdb.fastq.parse_fastq_filename])
        if labxdb.fastq.check_fastqs(fastqs) is False:
            raise Error('Name collision: Runs with the same name detected in different folder')

    irun = 1
    for name, path_list in fastqs.items():
        fastq_infos = labxdb.fastq.get_illumina_fastq_info(os.path.join(path_list[0]['path'], path_list[0]['fname']))
        # Max read length & Pair & Spots
        max_read_length = 0
        max_pair = 0
        spots = 0
        end = None
        for p in path_list:
            r = labxdb.fastq.get_illumina_fastq_info(os.path.join(p['path'], p['fname']), get_spots=True)
            if r['read_length'] > max_read_length:
                max_read_length = r['read_length']
            if r['pair'] > max_pair:
                max_pair = r['pair']
            # Keep the name of first end
            if end is None:
                end = p['end']
            # Add total spots for first end
            if p['end'] == end:
                spots += r['spots']
        # Paired
        if max_pair == 2:
            paired = True
        else:
            paired = False
        if dry_run:
            logger.info(f"DRY RUN: Adding {name} ({fastq_infos['flowcell']},{fastq_infos['barcode']},{max_read_length},paired:{paired})")
        else:
            logger.info(f'Adding run {name}')
            dbl.post('run/new', json=[{'run_ref':f'{ref_prefix}{irun:03}', 'run_order':1, 'tube_label':name, 'barcode':fastq_infos['barcode'], 'failed':False, 'flowcell':fastq_infos['flowcell'], 'paired':paired, 'max_read_length':max_read_length, 'spots':spots}])
        irun += 1

def import_raw_read_files(bulk, path_seq_raw, with_second_barcode=False, input_run_refs=[], exclude_run_refs=[], path_seq_run='.', print_summary=True, fastq_exts=['.fastq'], no_readonly=False, dry_run=False, dbl=None, config=None, num_processor=1, logger=None):
    # Parameters
    if logger is None:
        import logging as logger

    logger.info('Importing runs')
    fastqs = labxdb.fastq.find_fastqs(os.path.join(path_seq_raw, bulk), fastq_exts, [labxdb.fastq.parse_illumina_fastq_filename, labxdb.fastq.parse_fastq_filename])
    if labxdb.fastq.check_fastqs(fastqs) is False:
        raise Error('Name collision: Runs with the same name detected in different folder')

    # Create links to FASTQ files for each run
    done_runs = []
    created_folders = []
    for name, path_list in fastqs.items():
        # Get run info
        if dbl:
            search_criterion = []
            # Flowcell
            flowcells = set()
            for p in path_list:
                r = labxdb.fastq.get_illumina_fastq_info(os.path.join(p['path'], p['fname']))
                flowcells.add(r['flowcell'])
            if len(flowcells) == 1:
                flowcell = list(flowcells)[0]
                search_criterion.append('3 flowcell FUZZY '+flowcell)
            else:
                raise Error('Multiple flowcell detected in a single run')
            # Second barcode
            if with_second_barcode:
                name_parts = name.split('-')
                name = '-'.join(name_parts[:-1])
                second_barcode = name_parts[-1]
                search_criterion.append('3 second_barcode EQUAL '+second_barcode)
            # Tube label
            search_criterion.append('3 tube_label EQUAL '+name)
            # Query: Run
            runs = dbl.post('run', {'search_criterion':search_criterion, 'search_gate':'AND', 'limit':'ALL'})
            if len(input_run_refs) > 0:
                runs = [r for r in runs if r['run_ref'] in input_run_refs]
            if len(exclude_run_refs) > 0:
                runs = [r for r in runs if r['run_ref'] not in exclude_run_refs]
            if len(runs) == 0:
                if with_second_barcode:
                    logger.warning(f'{name}-{second_barcode} had no run')
                else:
                    logger.warning(f'{name} had no run')
                continue

        for r in runs:
            r['path_list'] = path_list

        # Create folders and symlinks
        for r in runs:
            full_path_run = os.path.join(path_seq_run, r['run_ref'])
            # Make folder
            if not os.path.exists(full_path_run):
                if dry_run:
                    logger.info(f'DRY RUN: os.makedir({full_path_run})')
                else:
                    os.mkdir(full_path_run)
                    created_folders.append(full_path_run)
            # Make symlinks
            for p in r['path_list']:
                fq_new_path = os.path.join(full_path_run, p['fname'])
                fq_rel_path = os.path.relpath(os.path.join(p['path'], p['fname']), full_path_run)
                if dry_run:
                    logger.info(f'DRY RUN: os.symlink({fq_rel_path}, {fq_new_path})')
                else:
                    if not os.path.exists(fq_new_path):
                        os.symlink(fq_rel_path, fq_new_path)
            # Done
            r['done'] = True

        done_runs.append(runs)

    # Change folder permissions to read-only
    if no_readonly is False:
        for f in created_folders:
            os.chmod(f, 0o0555)

    if print_summary:
        print('\nSummary')
        for runs in done_runs:
            for r in runs:
                if r['done']:
                    print(f"{r['tube_label']:<30}{r['flowcell']:<30}{r['barcode']:<30}{r['run_ref']:<15}{r['path_list'][0]['path']:<30}")
        print()

def check_exe(names):
    for name in names:
        if shutil.which(name) == None:
            raise Error(f'{name} missing')

def get_first_key(l, d):
    for k in l:
        if k in d:
            return d[k]

def main(argv=None):
    if argv is None:
        argv = sys.argv
    parser = argparse.ArgumentParser(description='Download and import sequencing read data.')
    # General
    group = parser.add_argument_group('Main')
    group.add_argument('-b', '--bulk', dest='bulk', action='store', help='Sequencing bulk name (defaults to "date_time")')
    group.add_argument('-y', '--format', dest='format', action='store', default='run_zstd', help='Organization/format of raw files. Supported: raw, run_zstd')
    group.add_argument('-r', '--path_seq_raw', dest='path_seq_raw', action='store', help='Path to raw seq.')
    group.add_argument('-l', '--path_seq_run', dest='path_seq_run', action='store', help='Path to run.')
    group.add_argument('-n', '--path_seq_tmp', dest='path_seq_tmp', action='store', help='Path to temporary.')
    group.add_argument('-k', '--no_readonly', dest='no_readonly', action='store_true', help='No read-only chmod')
    group.add_argument('-d', '--dry_run', dest='dry_run', action='store_true', help='Dry run')
    group.add_argument('-p', '--processor', dest='num_processor', action='store', type=int, default=6, help='Number of processor')
    group.add_argument('--path_config', dest='path_config', action='store', help='Path to config')
    group.add_argument('--fastq_exts', dest='fastq_exts', action='store', default='.fastq,.fastq.gz,.fastq.zst', help='FASTQ file extensions (comma separated).')
    group.add_argument('--http_url', '--labxdb_http_url', dest='labxdb_http_url', action='store', help='Database HTTP URL')
    group.add_argument('--http_login', '--labxdb_http_login', dest='labxdb_http_login', action='store', help='Database HTTP login')
    group.add_argument('--http_password', '--labxdb_http_password', dest='labxdb_http_password', action='store', help='Database HTTP password')
    group.add_argument('--http_path', '--labxdb_http_path', dest='labxdb_http_path', action='store', help='Database HTTP path')
    group.add_argument('--http_db', '--labxdb_http_db', dest='labxdb_http_db', action='store', help='Database HTTP DB')
    # Download
    group = parser.add_argument_group('Download')
    group.add_argument('-m', '--url_remote', dest='url_remote', action='append', help='Path to remote seq.')
    group.add_argument('--download_program', dest='download_program', action='store', default='wget', help='Tool used to download: wget')
    # Format
    group = parser.add_argument_group('Format')
    group.add_argument('--delete_download', dest='delete_download', action='store_true', help='Delete download folder.')
    group.add_argument('--squashfs_download', dest='squashfs_download', action='store_true', help='Create squashfs image of download folder.')
    group.add_argument('--flowcell_dir', dest='flowcell_dir', action='store_true', help='Use flowcell name as directory name.')
    # Import in DB
    group = parser.add_argument_group('DB import')
    group.add_argument('-e', '--ref_prefix', dest='ref_prefix', action='store', default='TMP_', help='Temporary run reference prefix.')
    # Import
    group = parser.add_argument_group('Import')
    group.add_argument('-t', '--input_run_refs', dest='input_run_refs', action='store', help='Only import these runs (comma separated references)')
    group.add_argument('--exclude_run_refs', dest='exclude_run_refs', action='store', help='Exclude these runs from import (comma separated references)')
    group.add_argument('--with_second_barcode', dest='with_second_barcode', action='store_true', help='Use the second barcode to identify FASTQ file(s).')
    # Step
    group = parser.add_argument_group('Step')
    group.add_argument('-g', '--make_download', dest='make_download', action='store_true', help='Download')
    group.add_argument('-z', '--make_format', dest='make_format', action='store_true', help='Format')
    group.add_argument('-s', '--make_staging', dest='make_staging', action='store_true', help='Input runs in staging tables.')
    group.add_argument('-i', '--make_import', dest='make_import', action='store_true', help='Import')
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
        if a in ['fastq_exts', 'input_run_refs', 'exclude_run_refs']:
            if v is None:
                config[a] = []
            else:
                config[a] = [r.strip() for r in v.split(',')]

    # Bulk option
    if config['make_download'] and 'bulk' not in config:
        date = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        if os.path.exists(os.path.join(get_first_key(['path_seq_tmp', 'path_seq_raw'], config), date)):
            print(f'ERROR: {date} already exists')
            return 1
        else:
            config['bulk'] = date
    elif 'bulk' not in config:
        print('ERROR: --bulk is required')
        return 1

    # Checks
    if config['make_download'] or config['make_format'] or config['make_staging'] or config['make_import']:
        # Check option
        if 'path_seq_raw' not in config:
            print('ERROR: --path_seq_raw is required')
            return 1
        # Check folders
        if not os.path.exists(config['path_seq_raw']):
            print(f"ERROR: {config['path_seq_raw']} not found")
            return 1
        if 'path_seq_tmp' in config and not os.path.exists(config['path_seq_tmp']):
            print(f"ERROR: {config['path_seq_tmp']} not found")
            return 1
    if config['make_import']:
        # Check option
        if 'path_seq_run' not in config:
            print('ERROR: --path_seq_run is required')
            return 1
        # Check folder
        if not os.path.exists(config['path_seq_run']):
            print(f"ERROR: {config['path_seq_run']} not found")
            return 1
    if config['make_download']:
        # Check for download software
        check_exe([config['download_program']])
    if config['make_format']:
        # Check for format software
        if config['format'].find('zstd') != -1:
            check_exe(['zstd'])
        if config['squashfs_download']:
            check_exe(['mksquashfs'])

    # Format
    path_format = os.path.join(get_first_key(['path_seq_tmp', 'path_seq_raw'], config), config['bulk'])
    path_download = os.path.join(path_format, 'download')
    if not config['dry_run'] and (config['make_download'] or config['make_format']):
        for p in [path_format, path_download]:
            if not os.path.exists(p):
                os.mkdir(p)
        log_filename = os.path.join(path_download, 'format.log')
    else:
        log_filename = None

    # Check permissions
    if not config['make_download'] and config['make_format'] and not os.access(path_download, os.W_OK):
        print(f"ERROR: {path_download} not writable")
        return 1

    # Logging
    logger = pfu.log.define_root_logger(f"load_{config['bulk']}", level='info', filename=log_filename)
    logger.info('Starting')

    try:
        # Download
        if config['make_download']:
            if 'url_remote' not in config:
                raise Error('Missing remote path')
            download_read_files(config['url_remote'], path_download, config, config['make_format'], download_program=config['download_program'], no_readonly=config['no_readonly'], dry_run=config['dry_run'], logger=logger)

        # Import
        fastqs = None
        if config['make_staging'] or config['make_import']:
            if 'labxdb_http_path' not in config and 'labxdb_http_db' not in config:
                if 'labxdb_http_path_seq' in config:
                    config['labxdb_http_path'] = config['labxdb_http_path_seq']
                else:
                    config['labxdb_http_db'] = 'seq'
            dbl = labxdb.DBLink(config.get('labxdb_http_url'), config.get('labxdb_http_login'), config.get('labxdb_http_password'), config.get('labxdb_http_path'), config.get('labxdb_http_db'))
        if config['make_format']:
            fastqs = format_raw_read_files(config['bulk'], config.get('path_seq_tmp'), config['path_seq_raw'], do_format=config['format'], delete_download=config['delete_download'], squashfs_download=config['squashfs_download'], flowcell_dir=config['flowcell_dir'], fastq_exts=config['fastq_exts'], no_readonly=config['no_readonly'], dry_run=config['dry_run'], num_processor=config['num_processor'], logger=logger)
        if config['make_staging']:
            import_staging(config['bulk'], config['path_seq_raw'], fastqs=fastqs, ref_prefix=config['ref_prefix'], fastq_exts=config['fastq_exts'], dry_run=config['dry_run'], dbl=dbl, logger=logger)
        if config['make_import']:
            import_raw_read_files(config['bulk'], config['path_seq_raw'], config['with_second_barcode'], config['input_run_refs'], config['exclude_run_refs'], config['path_seq_run'], fastq_exts=config['fastq_exts'], no_readonly=config['no_readonly'], dry_run=config['dry_run'], dbl=dbl, config=config, num_processor=config['num_processor'], logger=logger)
    except Error as e:
        logger.error(e.message)
        return 1

if __name__ == '__main__':
    sys.exit(main())
