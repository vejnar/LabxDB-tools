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
import csv
import json
import os
import sys

import labxdb

def uniquify(l):
    if len(l) == 0:
        return []
    else:
        u = [l[0]]
        for i, e in enumerate(l[1:]):
            if e != l[i]:
                u.append(e)
        return u

def xstr(s):
    if s is None:
        return ''
    return str(s)

# -----------------------
# Sample functions

def sample_name(infos, filters):
    n = f"{infos['project']['label_short']} - {infos['replicate']['label_short']} {infos['replicate']['replicate_ref']}"
    for k, v in filters:
        if n.find(k) != -1:
            n = n.replace(k, v)
            break
    return n

def sample_short_name(infos, filters):
    n = f"{infos['replicate']['label_short']} {infos['replicate']['replicate_ref']}"
    for k, v in filters:
        if n.find(k) != -1:
            n = n.replace(k, v)
            break
    return n

def sample_title(infos, filters):
    n = f"{infos['project']['label_short']} - {infos['replicate']['label_long']}"
    for k, v in filters:
        if n.find(k) != -1:
            n = n.replace(k, v)
            break
    return n

def project_label_long(infos, filters):
    n = f"{infos['project']['label_long']}"
    for k, v in filters:
        if n.find(k) != -1:
            n = n.replace(k, v)
            break
    return n

def project_label_short(infos, filters):
    n = f"{infos['project']['label_short']}"
    for k, v in filters:
        if n.find(k) != -1:
            n = n.replace(k, v)
            break
    return n

def sample_short_title(infos, filters):
    n = infos['replicate']['label_long']
    for k, v in filters:
        if n.find(k) != -1:
            n = n.replace(k, v)
            break
    return n

def species(infos):
    if infos['sample']['species'] == 'danRer':
        return 'Danio rerio'
    else:
        return infos['sample']['species']

def age(infos):
    if infos['sample']['age_hpf'] is not None:
        return float(infos['sample']['age_hpf'])
    else:
        return 'not applicable'

def tissue(infos):
    if infos['sample']['tissue'] is None:
        return 'embryo'
    else:
        return infos['sample']['tissue']

def genotype(infos):
    if infos['sample']['genotype'] != 'WT':
        return infos['sample']['genotype']
    else:
        return None

def second_barcode(infos):
    barcodes = list(set([run['second_barcode'] for run in infos['runs'] if run['second_barcode'] is not None]))
    if len(barcodes) > 0 :
        return ''.join(barcodes)
    else:
        return None

def get_sample_fields(label_filters, multiplex):
    c = {'sample_name': {'fn': sample_name, 'arg': label_filters},
         'sample_title': {'fn': sample_title, 'arg': label_filters},
         'bioproject_accession': {'default': None},
         'organism': {'fn': species},
         'strain': {'default': 'TU/AB'},
         'isolate': {'default': None},
         'breed': {'default': None},
         'cultivar': {'default': None},
         'ecotype': {'default': None},
         'age': {'fn': age},
         'dev_stage': {'column': ('sample', 'stage')},
         'sex': {'default': 'pooled male and female'},
         'tissue': {'fn': tissue},
         'biomaterial_provider': {'default': None},
         'birth_date': {'default': None},
         'birth_location': {'default': None},
         'breeding_history': {'default': None},
         'breeding_method': {'default': None},
         'cell_line': {'default': None},
         'cell_subtype': {'default': None},
         'cell_type': {'default': None},
         'collected_by': {'default': None},
         'collection_date': {'default': None},
         'culture_collection': {'default': None},
         'death_date': {'default': None},
         'disease': {'default': None},
         'disease_stage': {'default': None},
         'genotype': {'fn': genotype},
         'strain_maternal': {'column': ('sample', 'strain_maternal')},
         'strain_paternal': {'column': ('sample', 'strain_paternal')},
         'geo_loc_name': {'default': None},
         'growth_protocol': {'default': None},
         'health_state': {'default': None},
         'isolation_source': {'default': None},
         'lat_lon': {'default': None},
         'phenotype': {'default': None},
         'sample_type': {'default': None},
         'specimen_voucher': {'default': None},
         'store_cond': {'default': None},
         'stud_book_number': {'default': None},
         'treatment': {'column': ('sample', 'treatment')},
         'description': {'default': None},
         'molecule': {'column': ('sample', 'molecule')},
         'selection': {'column': ('sample', 'selection')},
         'condition': {'column': ('sample', 'condition')},
         'sample_ref': {'column': ('sample', 'sample_ref')},
         'replicate_ref': {'column': ('replicate', 'replicate_ref')},
         'replicate_order': {'column': ('replicate', 'replicate_order')},
         'project_label_long': {'fn': project_label_long, 'arg': label_filters},
         'project_label_short': {'fn': project_label_short, 'arg': label_filters},
         'sample_label_short': {'column': ('sample', 'label_short')},
         'replicate_label_short': {'column': ('replicate', 'label_short')}}
    if multiplex:
        c['barcode'] = {'fn': second_barcode}
    return c

def get_multiplex_sample_fields(label_filters, multiplex):
    c = get_sample_fields(label_filters, multiplex)
    c['sample_name'] = {'fn': sample_short_name, 'arg': label_filters}
    c['sample_title'] = {'fn': sample_short_title, 'arg': label_filters}
    return c

# -----------------------
# Data functions

def library_strategy(infos):
    if infos['sample']['library_protocol'] == 'dUTP':
        return 'RNA-seq'
    else:
        return 'OTHER'

def library_source(infos):
    if infos['sample']['molecule'] == 'DNA':
        return 'GENOMIC'
    elif infos['sample']['molecule'].find('RNA') != -1:
        return 'TRANSCRIPTOMIC'
    else:
        return 'OTHER'

def library_layout(infos):
    if infos['run']['paired']:
        return 'paired'
    else:
        return 'single'

def platform(infos):
    if infos['run']['platform'].find('Illumina') != -1:
        return 'ILLUMINA'

def filename(infos):
    return infos['run']['run_ref']+'_R1.fastq.gz'

def filename2(infos):
    if infos['run']['paired']:
        return infos['run']['run_ref']+'_R2.fastq.gz'
    else:
        return None

def get_data_fields(label_filters):
    return {'sample_name': {'fn': sample_name, 'arg': label_filters},
            'bioproject_accession': {'default': None},
            'biosample_accession': {'default': None},
            'title': {'fn': sample_title, 'arg': label_filters},
            'library_ID': {'column': ('run', 'run_ref')},
            'design_description': {'column': ('sample', 'molecule')},
            'library_strategy': {'fn': library_strategy},
            'library_source': {'fn': library_source},
            'library_selection': {'default': 'unspecified'},
            'library_layout': {'fn': library_layout},
            'platform': {'fn': platform},
            'instrument_model': {'column': ('run', 'platform')},
            'filetype': {'default': 'fastq'},
            'filename': {'fn': filename},
            'filename2': {'fn': filename2}}

def get_multiplex_data_fields(label_filters):
    c = get_data_fields(label_filters)
    c['sample_name'] = {'fn': sample_short_name, 'arg': label_filters}
    c['title'] = {'fn': sample_short_title, 'arg': label_filters}
    return c

def get_record(fields, infos):
    record = []
    for field_name, field in fields.items():
        val = None
        if 'default' in field:
            val = field['default']
        elif 'column' in field:
            val = infos[field['column'][0]][field['column'][1]]
        elif 'fn' in field:
            if 'arg' in field:
                val = field['fn'](infos, field['arg'])
            else:
                val = field['fn'](infos)
        record.append(val)
    return record

def get_multiplex_record(fields, multiplex_samples, flowcell_alias, tube_label, unexported_barcodes):
    combined_record = []
    for ifield, field_name in enumerate(fields.keys()):
        mfield = [m['record'][ifield] for m in multiplex_samples]
        uniq = uniquify(mfield)
        if field_name == 'sample_name' or field_name == 'title' or field_name == 'sample_title':
            new_name = 'Raw multiplex: ' + ';'.join(mfield + ['unrelated']*len(unexported_barcodes))
        elif (field_name == 'filename' or field_name == 'filename2') and any(uniq):
            if field_name == 'filename':
                new_name = f'{flowcell_alias}_{tube_label}_R1.fastq.gz'
            elif field_name == 'filename2':
                new_name = f'{flowcell_alias}_{tube_label}_R2.fastq.gz'
        elif field_name == 'library_ID':
            new_name = ';'.join([xstr(u) for u in uniq + ['unrelated']*len(unexported_barcodes)])
        elif field_name == 'barcode':
            new_name = ';'.join([xstr(u) for u in uniq + unexported_barcodes])
        else:
            new_name = ';'.join([xstr(u) for u in uniq])
        combined_record.append(new_name)
    return combined_record

def merge_refs(infos, column):
    return ','.join([run[column] for run in infos['runs'] if run[column] is not None])

def get_exported_fields(label_filters):
    return {'project_ref': {'column': ('project', 'project_ref')},
            'sample_ref': {'column': ('sample', 'sample_ref')},
            'replicate_ref': {'column': ('replicate', 'replicate_ref')},
            'replicate_order': {'column': ('replicate', 'replicate_order')},
            'project_label_short': {'column': ('project', 'label_short')},
            'label_short': {'column': ('replicate', 'label_short')},
            'label_long': {'column': ('replicate', 'label_long')},
            'age': {'fn': age},
            'dev_stage': {'column': ('sample', 'stage')},
            'treatment': {'column': ('sample', 'treatment')},
            'selection': {'column': ('sample', 'selection')},
            'condition': {'column': ('sample', 'condition')},
            'publication_ref': {'column': ('replicate', 'publication_ref')},
            'replicate_sra_ref': {'column': ('replicate', 'sra_ref')},
            'run_refs': {'fn': merge_refs, 'arg': 'run_ref'},
            'run_sra_refs': {'fn': merge_refs, 'arg': 'sra_ref'}}

# -----------------------
# Main

def main(argv=None):
    if argv is None:
        argv = sys.argv
    parser = argparse.ArgumentParser(description='Export data to SRA.')
    parser.add_argument('-r', '--replicates', dest='replicates', action='store', help='Replicates (comma separated)')
    parser.add_argument('--path_replicates', dest='path_replicates', action='store', help='Path to replicates')
    parser.add_argument('-f', '--path_label_filters', dest='path_label_filters', action='store', help='Path to label filters')
    parser.add_argument('-m', '--export_multiplex', dest='export_multiplex', action='store_true', default=False, help='Export multiplex data')
    parser.add_argument('-e', '--exclude_exported', dest='exclude_exported', action='store_true', default=False, help='Exclude already exported data')
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
        if a in ['replicates']:
            if v is None:
                config[a] = []
            else:
                config[a] = [r.strip() for r in v.split(',')]

    # Init. DBLink
    if 'labxdb_http_path' not in config and 'labxdb_http_db' not in config:
        if 'labxdb_http_path_seq' in config:
            config['labxdb_http_path'] = config['labxdb_http_path_seq']
        else:
            config['labxdb_http_db'] = 'seq'
    dbl = labxdb.DBLink(config.get('labxdb_http_url'), config.get('labxdb_http_login'), config.get('labxdb_http_password'), config.get('labxdb_http_path'), config.get('labxdb_http_db'))

    # Init. CSV export
    csv.register_dialect('excel-tab-unix', delimiter='\t', lineterminator='\n', quoting=csv.QUOTE_MINIMAL)

    # Replicates
    replicates = config['replicates']
    if 'path_replicates' in config:
        with open(config['path_replicates']) as f:
            replicates.extend([l.split(',')[0] for l in f.read().strip().split()])

    # Check replicates
    if len(replicates) == 0:
        print('ERROR: No replicate defined')
        return 1

    # Get infos
    project_infos = {}
    sample_infos = {}
    replicate_infos = []
    for replicate_ref in replicates:
        # Query
        replicate = dbl.get('replicate/get-ref/'+replicate_ref)[0][0]
        sample = dbl.get('sample/get-ref/'+replicate['sample_ref'])[0][0]
        project = dbl.get('project/get-ref/'+sample['project_ref'])[0][0]
        # Add project ref to replicate
        replicate['project_ref'] = sample['project_ref']
        # Save
        replicate_infos.append(replicate)
        sample_infos[sample['sample_ref']] = sample
        project_infos[project['project_ref']] = project
    replicate_infos.sort(key=lambda r: (r['project_ref'], r['sample_ref'], r['replicate_order']))

    # Label filters
    if 'path_label_filters' in config:
        with open(config['path_label_filters']) as fcsv:
            label_filters = [row for row in csv.reader(fcsv)]
    else:
        label_filters = []

    # SRA fields
    sample_fields = get_sample_fields(label_filters, config['export_multiplex'])
    data_fields = get_data_fields(label_filters)
    exported_fields = get_exported_fields(label_filters)
    # SRA fields for multiplex
    if config['export_multiplex']:
        multiplex_sample_fields = get_multiplex_sample_fields(label_filters, config['export_multiplex'])
        multiplex_data_fields = get_multiplex_data_fields(label_filters)

    # Get records
    sample_records = []
    data_records = []
    exported_records = []
    data_files = []
    multiplex_samples_records = {}
    multiplex_runs = {}
    for replicate in replicate_infos:
        # Get run(s)
        runs = dbl.post('run', {'search_criterion':['3 failed EQUAL FALSE', '3 replicate_ref EQUAL '+replicate['replicate_ref']], 'search_gate':'AND', 'limit':'ALL'})

        # Already exported run
        if config['exclude_exported'] and replicate['sra_ref'] is not None:
            print(f"WARNING: {replicate['replicate_ref']} \"{replicate['label_short']}\" already exported in {replicate['publication_ref']}")
            erecord = get_record(exported_fields, {'runs':runs, 'replicate':replicate, 'sample':sample_infos[replicate['sample_ref']], 'project':project_infos[replicate['project_ref']]})
            exported_records.append(erecord)
            continue

        # SRA sample
        srecord = get_record(sample_fields, {'runs':runs, 'replicate':replicate, 'sample':sample_infos[replicate['sample_ref']], 'project':project_infos[replicate['project_ref']]})
        sample_records.append(srecord)

        # Multiplex sample
        if config['export_multiplex']:
            multiplex_samples_records[replicate['replicate_ref']] = get_record(multiplex_sample_fields, {'runs':runs, 'replicate':replicate, 'sample':sample_infos[replicate['sample_ref']], 'project':project_infos[replicate['project_ref']]})

        # SRA data
        for run in runs:
            # Data record
            drecord = get_record(data_fields, {'run':run, 'replicate':replicate, 'sample':sample_infos[replicate['sample_ref']], 'project':project_infos[replicate['project_ref']]})
            data_records.append(drecord)
            # Data
            data_files.append({'type':'single', 'run_ref':run['run_ref'], 'pattern':'_R1.fastq'})
            if run['paired']:
                data_files.append({'type':'single', 'run_ref':run['run_ref'], 'pattern':'_R2.fastq'})

            # Multiplex run
            if config['export_multiplex'] and run['second_barcode'] is not None:
                key = (run['flowcell'], run['tube_label'])
                mdrecord = get_record(multiplex_data_fields, {'run':run, 'replicate':replicate, 'sample':sample_infos[replicate['sample_ref']], 'project':project_infos[replicate['project_ref']]})
                if key in multiplex_runs:
                    multiplex_runs[key].append({'record':mdrecord, 'run':run})
                else:
                    multiplex_runs[key] = [{'record':mdrecord, 'run':run}]

    # Append multiplex samples
    if config['export_multiplex']:
        multiplex_samples = {}
        for k, mruns in multiplex_runs.items():
            flowcell, tube_label = k
            flowcell_alias = flowcell.split(':')[-1]

            # Get unexported barcodes
            runs = dbl.post('run', {'search_criterion':['3 failed EQUAL FALSE', '3 flowcell EQUAL '+flowcell, '3 tube_label EQUAL '+tube_label], 'search_gate':'AND', 'limit':'ALL'})
            barcodes = set([r['second_barcode'] for r in runs])
            exported_barcodes = [m['run']['second_barcode'] for m in mruns]
            unexported_barcodes = list(barcodes.difference(set(exported_barcodes)))

            # Data record
            data_records.append(get_multiplex_record(data_fields, mruns, flowcell_alias, tube_label, unexported_barcodes))
            # Data
            data_files.append({'type':'multiplex', 'flowcell':flowcell_alias, 'tube_label':tube_label, 'pattern':'_R1.fastq'})
            if run['paired']:
                data_files.append({'type':'multiplex', 'flowcell':flowcell_alias, 'tube_label':tube_label, 'pattern':'_R2.fastq'})
            
            # Multiplex sample
            key = tuple([m['run']['replicate_ref'] for m in mruns])
            multiplex_samples[key] = {'unexported_barcodes':unexported_barcodes}

        for i, (msamples, minfos) in enumerate(multiplex_samples.items()):
            sample_records.append(get_multiplex_record(multiplex_sample_fields, [{'record':multiplex_samples_records[n]} for n in msamples], '', '', minfos['unexported_barcodes']))

    # Export
    for records, fields, name in [(sample_records, sample_fields, 'sra_samples.tsv'), (data_records, data_fields, 'sra_data.tsv'), (exported_records, exported_fields, 'sra_exported.tsv')]:
        with open(name, 'w', newline='') as f:
            writer = csv.writer(f, dialect='excel-tab-unix')
            writer.writerow(fields.keys())
            writer.writerows(records)
    with open('data.json', 'wt') as f:
        json.dump(data_files, f)

if __name__ == '__main__':
    sys.exit(main())
