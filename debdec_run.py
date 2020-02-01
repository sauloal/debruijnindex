#!/usr/bin/env python3

import sys
import os

from debdec import sliding_windows
from debdec_create import get_db

def readfile(k, fhd):
    p = True

    d = get_db(k)
    print(d)

    table     = {base:idx for idx, base in enumerate("ACGT")}
    table_inv = [c for c, i in sorted(table.items(), key=lambda x:x[1])]

    for rownum, line in enumerate(fhd):
        i     = [table[c] for c in line.rstrip()]
        kmers = sliding_windows(i, k)
        print(f"{' '*18}{line.strip()}")
        for knum, kmer in enumerate(kmers):
            v   = d(kmer)
            dna = "".join([table_inv[i] for i in kmer]) if p else ""
            print(f"{rownum:5d} {knum:5d} {v:5d} {' '*knum}{dna}")

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("usage: {} <k-mer length> [<filename>]".format(sys.argv[0]))
        exit(1)
    
    try:
        k = int(sys.argv[1])
    except:
        print("usage: {} <k-mer length> [<filename>]".format(sys.argv[0]))
        exit(1)
    
    if k < 1:
        print("k-mer length must be at least 1.")

    fhd = sys.stdin

    if len(sys.argv) == 3:
        filename = sys.argv[2]
        assert os.path.exists(filename)
        fhd = open(filename, 'r')

    readfile(k, fhd)

    if len(sys.argv) == 3:
        fhd.close()
