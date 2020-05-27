# -*- coding: utf-8 -*-

#
# Copyright (C) 2018-2020 Charles E. Vejnar
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
#

import os
import tempfile
import time
import urllib.request
import xml.etree.ElementTree as ET

def get_samples_infos(srp, url_search='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&retmax=1000&term=', url_get='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=', import_runs=[], save_sra_xml=False, verbose=False):
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
        sra_root = ET.parse(xmlo).getroot()
        # Check study
        study = sra_root.find('.//STUDY')
        if (srp == study.attrib['accession'] or srp == study.attrib['alias']) == False:
            time.sleep(1)
            continue
        # Sample
        sample = {'runs': []}
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
            paired = False
        elif sra_root.find('.//DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT/PAIRED') is not None:
            paired = True
        else:
            paired = None
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
                   'paired': paired,
                   'platform': platform}
            sra_stat_per_read = {c.get('index'): c.get('count') for c in sra_run.find('./Statistics')}
            if '0' in sra_stat_per_read:
                run['nread1'] = sra_stat_per_read['0']
            if '1' in sra_stat_per_read:
                run['nread2'] = sra_stat_per_read['1']
            if len(import_runs) == 0 or (len(import_runs) > 0 and run['ref'] in import_runs):
                sample['runs'].append(run)
        if len(sample['runs']) > 0:
            sample['runs'].sort(key=lambda x: x['ref'])
            samples.append(sample)
        # Sleep to reduce request rate
        if i % 3 == 0 or i % 10 == 0:
            time.sleep(1)
    samples.sort(key=lambda x: x['ref'])
    return srp_title, samples
