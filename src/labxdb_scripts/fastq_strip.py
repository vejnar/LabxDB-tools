#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# Copyright (C) 2017-2022 Charles E. Vejnar
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
#

import sys

def main(argv=None):
    if argv is None:
        argv = sys.argv
    # Arguments
    if len(argv) != 3:
        print('ERROR: missing argument')
        return 1
    else:
        pin = argv[1]
        pout = argv[2]

    with open(pin, 'rt', buffering=4096000) as fin, open(pout, 'wt', buffering=4096000) as fout:
        nrecord = 0
        mlen = 0
        record = []
        for l in fin:
            record.append(l)
            if len(record) == 4:
                fout.write(record[0].split(' ')[0] + '\n')
                fout.write(record[1])
                fout.write('+\n')
                fout.write(record[3])
                if len(record[1].strip()) > mlen:
                    mlen = len(record[1].strip())
                record = []
                nrecord += 1
        print(nrecord, mlen)

if __name__ == '__main__':
    sys.exit(main())
