# -*- coding: utf-8 -*-

#
# Copyright © 2018 Charles E. Vejnar
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

def parse_no_end_fastq_filename(fname):
    m = re.match('(.+)\.fastq', fname)
    if m:
        return {'name':m.group(1), 'end':None}
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

def get_illumina_fastq_info(line):
    """Queries the first entry in the fastq file and reports information, assuming Illumina semantics:

    instrument, run number, flowcell, lane, pair number, read length, barcode
    @DCM97JN1:188:C11R8ACXX:5:1101:1208:2164 1:N:0:GCTCATGA
    """
    try:
        info = {}

        # Get info from ID
        seq_ids = line.split()
        if len(seq_ids) == 2:
            identifier, meta = seq_ids
        else:
            identifier, meta = seq_ids[0], None
        id_fields = identifier.split(':')

        info['instrument'] = id_fields[0][1:]
        info['run number'] = id_fields[1]
        info['flowcell'] = id_fields[2]
        info['lane'] = int(id_fields[3])
        if meta:
            info['pair'] = int(meta[0])
            info['barcode'] = meta.split(':')[-1]
        else:
            info['pair'] = int(id_fields[-1][0])
            info['barcode'] = None

        return info
    except:
        return None

def get_pacbio_fastq_info(line):
    """Queries the first entry in the fastq file and reports information, assuming PacBio semantics:

    {movieName}/{holeNumber}/ccs
    """

    try:
        info = {}

        # Get info from ID
        seq_ids = line.split('/')
        if len(seq_ids) >= 3:
            info['flowcell'] = seq_ids[0][1:]

        return info
    except:
        return None

def get_fastq_info(path_seq, fastq_fns=[get_illumina_fastq_info, get_pacbio_fastq_info], get_spots=False):
    # Get first read ID
    if path_seq.endswith('.gz'):
        fqf = gzip.open(path_seq, 'rt')
    elif path_seq.endswith('.zst'):
        fqf = pfu.zstd.open(path_seq, 'rt')
    else:
        fqf = open(path_seq, 'rt', buffering=1024*1024)
    line = fqf.readline().rstrip('\n\r')
    assert line.startswith('@'), f'Invalid FASTQ header: {line} in {path_seq}'

    # Info
    info = {}
    for fastq_fn in fastq_fns:
        m = fastq_fn(line)
        if m is not None:
            info = m
            break

    # Get first sequence
    seq = fqf.readline().rstrip('\n\r')
    info['read_length'] = len(seq)

    # Get number of spots
    if get_spots:
        nread = 1 # One read already
        max_len = info['read_length']
        while True:
            try:
                read = [next(fqf) for i in range(4)]
            except StopIteration:
                break
            l = len(read[1]) - 1
            if l > max_len:
                max_len = l
            nread += 1
        info['spots'] = nread
        info['read_length'] = max_len

    # Close FASTQ
    fqf.close()

    return info
