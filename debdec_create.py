#!/usr/bin/env python3

import sys

from debdec import make_decoder

def get_db(k):
    return make_decoder(k)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("usage: {} <k-mer length>".format(sys.argv[0]))
        exit(1)
    
    try:
        k = int(sys.argv[1])
    except:
        print("usage: {} <k-mer length>".format(sys.argv[0]))
        exit(1)
    
    if k < 1:
        print("k-mer length must be at least 1.")

    get_db(k)
