# -*- coding: utf-8 -*-

#
# Copyright (C) 2018-2022 Charles E. Vejnar
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
#

import json
import os
import tempfile
import time
import urllib.request
import xml.etree.ElementTree as ET

def get_samples_infos(srp, url_search='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&retmax=1000&term=', url_get='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=', import_runs=[], save_sra_xml=False, save_project_json=False, verbose=False):
    samples = []
    srp_title = None
    rep_search = urllib.request.urlopen(url_search+srp)
    for i, sra_id in enumerate(ET.parse(rep_search).getroot().iter('Id')):
        if verbose:
            print('Downloading infos #', i, 'ID', sra_id.text)
        rep_get = urllib.request.urlopen(url_get+sra_id.text)
        # Read XML
        if save_sra_xml:
            xmlo = os.path.join(tempfile.gettempdir(), '%s_%s.xml'%(srp, i))
            with open(xmlo, 'wb') as f:
                f.write(rep_get.read())
        else:
            xmlo = rep_get
        # Parse
        title, sample = parse_sra_xml(xmlo, srp, import_runs)
        if title is not None:
            srp_title = title
        if sample is not None and len(sample['runs']) > 0:
            samples.append(sample)
        # Sleep to reduce request rate
        time.sleep(1)
    samples.sort(key=lambda x: (x['ref'], x['runs'][0]['ref']))
    # Project
    project = {'title': srp_title, 'samples': samples}
    if save_project_json:
        with open(os.path.join(tempfile.gettempdir(), f'{srp}.json'), 'wt') as fout:
            json.dump(project, fout)
    return project

def parse_sra_xml(f, srp=None, import_runs=[]):
    sra_root = ET.parse(f).getroot()
    # Sample
    sample = {'runs': []}
    # Check study
    study = sra_root.find('.//STUDY')
    if srp is not None and (srp == study.attrib['accession'] or srp == study.attrib['alias']) == False:
        return sample
    # Ref
    sra_sample = sra_root.find('.//SAMPLE')
    if sra_sample is not None:
        sample['ref'] = sra_sample.attrib['accession']
        sample['name'] = sra_sample.attrib['alias']
    # Title
    sra_title = sra_root.find('.//SAMPLE/TITLE')
    if sra_title is not None:
        sample['label'] = sra_title.text
    elif 'name' in sample:
        sample['label'] = sample['name']
    # Layout
    if sra_root.find('.//DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT/SINGLE') is not None:
        sample['paired'] = False
    elif sra_root.find('.//DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT/PAIRED') is not None:
        sample['paired'] = True
    else:
        sample['paired'] = None
    # Platform
    sra_illumina = sra_root.find('.//PLATFORM/ILLUMINA/INSTRUMENT_MODEL')
    sra_pacbio = sra_root.find('.//PLATFORM/PACBIO_SMRT/INSTRUMENT_MODEL')
    if sra_illumina is not None:
        platform = sra_illumina.text
    elif sra_pacbio is not None:
        platform = sra_pacbio.text
    else:
        platform = None
    # SRP title
    sra_title = sra_root.find('.//STUDY_TITLE')
    if sra_title is not None:
        srp_title = sra_title.text
    else:
        srp_title = None
    # Attribute(s)
    sample['attributes'] = [[sra_att.find('./TAG').text, sra_att.find('./VALUE').text] for sra_att in sra_root.iter('SAMPLE_ATTRIBUTE')]
    # Run(s)
    for sra_run in sra_root.iter('RUN'):
        # Platform
        sra_illumina = sra_run.find('./PLATFORM/ILLUMINA/INSTRUMENT_MODEL')
        sra_pacbio = sra_root.find('./PLATFORM/PACBIO_SMRT/INSTRUMENT_MODEL')
        if sra_illumina is not None:
            platform = sra_illumina.text
        elif sra_pacbio is not None:
            platform = sra_pacbio.text
        # Stats
        run = {'ref': sra_run.attrib['accession'],
               'spots': int(sra_run.attrib['total_spots']),
               'paired': False,
               'platform': platform}
        sra_stat_per_read = {c.get('index'): c.get('count') for c in sra_run.find('./Statistics')}
        if '0' in sra_stat_per_read:
            run['nread1'] = int(sra_stat_per_read['0'])
        if '1' in sra_stat_per_read:
            run['nread2'] = int(sra_stat_per_read['1'])
            run['paired'] = True
        # URL
        for sra_file in sra_run.iter('SRAFile'):
            if sra_file.attrib['cluster'] == 'public' and sra_file.attrib['supertype'] == 'Primary ETL':
                run['sra_url'] = sra_file.attrib['url']
                run['sra_size'] = int(sra_file.attrib['size'])
                md5 = sra_file.get('md5')
                if md5 is not None:
                    run['sra_md5'] = md5
        # Append run
        if len(import_runs) == 0 or (len(import_runs) > 0 and run['ref'] in import_runs):
            sample['runs'].append(run)
    if len(sample['runs']) > 0:
        sample['runs'].sort(key=lambda x: x['ref'])
    return srp_title, sample
