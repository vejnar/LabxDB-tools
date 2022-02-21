#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# Copyright (C) 2019-2022 Charles E. Vejnar
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
#

import argparse
import json
import os
import re
import sys
import tempfile

import labxdb
import labxdb.ncbi

def xstr(s):
    if s is None:
        return ''
    return str(s)

def search_replicate_ref(sra_samples):
    nwarn = 0
    for sra_sample in sra_samples:
        found = False
        name = ''
        for field in ['name', 'label']:
            if field in sra_sample:
                name = sra_sample[field]
        for k, v in sra_sample['attributes']:
            if k == 'replicate_ref':
                sra_sample['replicate_ref'] = v
                found = True
        if found is False:
            print('WARNING: no ID found in %s %s'%(name, [r['spots'] for r in sra_sample['runs']]))
            nwarn += 1
    return nwarn

def update_ref(sra_samples, publication_ref, dbl):
    queries = []
    nwarn = 0
    for sra_sample in sra_samples:
        if sra_sample['name'].startswith('Raw multiplex:'):
            print('Skipping multiplex', sra_sample['name'])
        elif 'replicate_ref' in sra_sample:
            # Get replicate
            result = dbl.get('replicate/get-ref/'+sra_sample['replicate_ref'])[0]
            if len(result) == 0:
                print('WARNING: replicate %s not found in DB'%sra_sample['replicate_ref'])
                nwarn += 1
                continue
            replicate = result[0]
            print('  Replicate\n  > Local: {: <10} {: <30} {: <20}'.format(replicate['replicate_ref'], replicate['label_short'], xstr(replicate['sra_ref'])))
            print('  > SRA:   {: <41} {: <20}'.format(sra_sample['name'], sra_sample['ref']))

            # SRA ref.
            if replicate['sra_ref'] is None:
                queries.append(['replicate/edit/%s'%replicate['replicate_id'], {'sra_ref':sra_sample['ref']}])
            else:
                print('WARNING: %s already attached to %s (%s) in DB'%(replicate['replicate_ref'], replicate['sra_ref'], replicate['publication_ref']))
                nwarn += 1

            # Publication ref.
            if replicate['publication_ref'] is None:
                queries.append(['replicate/edit/%s'%replicate['replicate_id'], {'publication_ref':publication_ref}])

            # Runs
            runs = dbl.post('run', {'search_criterion':['3 replicate_ref EQUAL '+replicate['replicate_ref']], 'limit':'ALL'})
            for sra_run in sra_sample['runs']:
                found = False
                already = False
                for run in runs:
                    if run['sra_ref'] is None:
                        if sra_run['spots'] == run['spots']:
                            found = True
                            break
                    else:
                        already = True
                if found:
                    print('    > Run: %s'%sra_run['ref'])
                    queries.append(['run/edit/%s'%run['run_id'], {'sra_ref':sra_run['ref']}])
                elif already:
                    print('WARNING: %s already attached to %s in DB'%(run['run_ref'], run['sra_ref']))
                    nwarn += 1
                else:
                    print('WARNING: no run found for %s'%sra_run['ref'])
                    nwarn += 1
    return queries, nwarn

def check_replicates(sra_samples, publication_ref, dbl):
    sra_runs = []
    nwarn = 0
    for sra_sample in sra_samples:
        replicates = dbl.post('replicate', {'search_criterion':['2 sra_ref FUZZY '+sra_sample['ref']], 'limit':'ALL'})
        if len(replicates) == 0:
            print('WARNING: %s not found in DB (%s)'%(sra_sample['ref'], sra_sample['label']))
            nwarn += 1
        elif len(replicates) > 1:
            print('WARNING: %s multiple replicates found in DB'%sra_sample['ref'])
            nwarn += 1
        else:
            replicate = replicates[0]
            if replicate['publication_ref'] == publication_ref:
                for sra_run in sra_sample['runs']:
                    sra_run['replicate_ref'] = replicate['replicate_ref']
                sra_runs.extend(sra_sample['runs'])
            else:
                print('WARNING: %s already attached to %s in DB'%(replicate['replicate_ref'], replicate['publication_ref']))
                nwarn += 1
    return sra_runs, nwarn

def check_runs(sra_runs, dbl):
    nwarn = 0
    nerror = 0
    for sra_run in sra_runs:
        sra_run['done'] = False
    for sra_run in sra_runs:
        if sra_run['done'] == False:
            runs = dbl.post('run', {'search_criterion':['3 replicate_ref EQUAL '+sra_run['replicate_ref'], '3 sra_ref FUZZY '+sra_run['ref']], 'limit':'ALL', 'search_gate':'AND'})
            if len(runs) == 0:
                print('WARNING: no run was found in DB for %s (%s)'%(sra_run['ref'], sra_run['spots']))
                nwarn += 1
            elif len(runs) > 1:
                print('WARNING: more than one run was found in DB for %s'%sra_run['ref'])
                nwarn += 1
            else:
                run = runs[0]
                # Get total spots in SRA
                sra_refs = run['sra_ref'].split(',')
                sra_spots = 0
                for sr in sra_runs:
                    if sr['ref'] in sra_refs and sr['done'] == False:
                        sra_spots += sr['spots']
                        sr['done'] = True
                print('Found run %s for %s'%(sra_run['ref'], run['run_ref']))
                if run['spots'] != sra_spots:
                    print('ERROR: %s %s != %s'%(sra_run['ref'], run['spots'], sra_spots))
                    nerror += 1
                    sys.exit()
    print('Missing run', [sra_run['ref'] for sra_run in sra_runs if sra_run['done'] == False])
    return nwarn, nerror

def main(argv=None):
    if argv is None:
        argv = sys.argv
    parser = argparse.ArgumentParser(description='Import SRA IDs to LabxDB seq.')
    parser.add_argument('-u', '--update', dest='update', action='store_true', help='Update.')
    parser.add_argument('-c', '--check', dest='check', action='store_true', help='Check.')
    parser.add_argument('-d', '--dry', dest='dry', action='store_true', help='Dry run.')
    parser.add_argument('-p', '--publication_ref', dest='publication_ref', action='store', required=True, help='Publication reference.')
    parser.add_argument('-x', '--save_sra_xml', dest='save_sra_xml', action='store_true', help='Save project XML from SRA.')
    parser.add_argument('-s', '--save_project_json', dest='save_project_json', action='store_true', help='Save project (JSON).')
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

    # Init. DBLink
    if 'labxdb_http_path' not in config and 'labxdb_http_db' not in config:
        if 'labxdb_http_path_seq' in config:
            config['labxdb_http_path'] = config['labxdb_http_path_seq']
        else:
            config['labxdb_http_db'] = 'seq'
    dbl = labxdb.DBLink(config.get('labxdb_http_url'), config.get('labxdb_http_login'), config.get('labxdb_http_password'), config.get('labxdb_http_path'), config.get('labxdb_http_db'))

    # Publication
    publication = dbl.get('publication/get-ref/'+config['publication_ref'])[0][0]
    print('Title:', publication['title'])
    
    # Add Publication to Option
    if not config['dry']:
        options = dbl.post('option', {'search_criterion':['0 group_name EQUAL publication_ref'], 'limit':'ALL'})
        if config['publication_ref'] not in [o['option'] for o in options]:
            print(f"Adding publication {config['publication_ref']} in Option")
            dbl.post('option/new', json=[{'group_name':'publication_ref', 'option':config['publication_ref']}])

    # Update & Check
    nwarn = 0
    nerror = 0
    for sra_ref in publication['sra_ref'].split(','):
        # Get infos
        print('Loading', sra_ref)
        path_info = os.path.join(tempfile.gettempdir(), sra_ref+'.json')
        if os.path.exists(path_info):
            project = json.load(open(path_info))
        else:
            project = labxdb.ncbi.get_samples_infos(sra_ref, save_sra_xml=config['save_sra_xml'], save_project_json=config['save_project_json'], verbose=True)
        # Update
        print('\n>', sra_ref)
        if config['update']:
            nwarn += search_replicate_ref(project['samples'])
            queries, nw = update_ref(project['samples'], config['publication_ref'], dbl)
            nwarn += nw
            for q in queries:
                if not config['dry']:
                    dbl.post(q[0], json=[q[1]])
        # Check
        if config['check']:
            sra_runs, nw = check_replicates(project['samples'], config['publication_ref'], dbl)
            nwarn += nw
            nw, ne = check_runs(sra_runs, dbl)
            nwarn += nw
            nerror += ne
    print(f'{nwarn} warning(s), {nerror} error(s)')

if __name__ == '__main__':
    sys.exit(main())
