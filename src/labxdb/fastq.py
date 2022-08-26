# -*- coding: utf-8 -*-

#
# Copyright (C) 2018-2022 Charles E. Vejnar
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
#

import gzip
import os
import re

import pyfnutils as pfu
import pyfnutils.zstd

def parse_illumina_fastq_filename(fname):
    """Illumina semantics: *_NUMBER_L00X_R1_YYY.fastq.XX"""
    m = re.match('(.+)_([a-zA-Z0-9]+)_L(\d{3})_([R,I][1,2,3])_\d{3}\.f', fname)
    if m:
        return {'name':m.group(1), 'barcode':m.group(2), 'lane':m.group(3), 'end':m.group(4)}
    else:
        return None

def parse_fastq_filename(fname):
    m = re.match('(.+)_([R,I][1,2,3])\.f', fname)
    if m:
        return {'name':m.group(1), 'end':m.group(2)}
    else:
        return None

def find_fastqs(path_seq, fastq_exts=['.fastq'], fastq_fns=[parse_fastq_filename]):
    fastqs = {}
    for path, dirs, files in os.walk(path_seq, followlinks=True):
        for fname in files:
            if any([fname.endswith(e) for e in fastq_exts]):
                for fastq_fn in fastq_fns:
                    m = fastq_fn(fname)
                    if m is not None:
                        m['path'] = path
                        m['fname'] = fname
                        if m['name'] not in fastqs:
                            fastqs[m['name']] = []
                        fastqs[m['name']].append(m)
                        break
    # Sort files
    for k, v in fastqs.items():
        v.sort(key=lambda x: x['fname'])
    # Sort keys
    sorted_fastqs = {}
    for k in sorted(fastqs.keys()):
        sorted_fastqs[k] = fastqs[k]
    return sorted_fastqs

def check_fastqs(fastqs):
    for name, path_list in fastqs.items():
        if len(set([p['path'] for p in path_list])) != 1:
            return False
    return True

def get_illumina_fastq_info(ilmn_fastq, get_spots=False):
    """Queries the first entry in the fastq file and reports information, assuming Illumina semantics:

    flowcell, lane, pair number, read length, barcode
    @DCM97JN1:188:C11R8ACXX:5:1101:1208:2164 1:N:0:GCTCATGA
    """
    # Get first read ID
    if ilmn_fastq.endswith('.gz'):
        fqf = gzip.open(ilmn_fastq, 'rt')
    elif ilmn_fastq.endswith('.zst'):
        fqf = pfu.zstd.open(ilmn_fastq, 'rt')
    else:
        fqf = open(ilmn_fastq, 'rt', buffering=1024*1024)
    line = fqf.readline().rstrip('\n\r')
    assert line.startswith('@'), f'Invalid FASTQ header: {line} in {ilmn_fastq}'
    seq_ids = line.split()

    # Info
    infos = {}

    # Get first sequence
    seq = fqf.readline().rstrip('\n\r')
    infos['read_length'] = len(seq)

    # Get info from ID
    if len(seq_ids) == 2:
        identifier, meta = seq_ids
    else:
        identifier, meta = seq_ids[0], None
    id_fields = identifier.split(':')
    infos['flowcell'] = ':'.join(id_fields[:-4])[1:]
    infos['lane'] = int(id_fields[-4])
    if meta:
        infos['pair'] = int(meta[0])
        infos['barcode'] = meta.split(':')[-1]
    else:
        infos['pair'] = int(id_fields[-1][0])
        infos['barcode'] = None

    # Get number of spots
    if get_spots:
        nread = 1 # One read already
        max_len = infos['read_length']
        while True:
            try:
                read = [next(fqf) for i in range(4)]
            except StopIteration:
                break
            l = len(read[1]) - 1
            if l > max_len:
                max_len = l
            nread += 1
        infos['spots'] = nread
        infos['read_length'] = max_len

    # Close FASTQ
    fqf.close()

    return infos
