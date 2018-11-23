#!/usr/bin/python3

import argparse
import math
import sys

parser = argparse.ArgumentParser(
    description='Create histogram bins for values in GNUplot data file')
parser.add_argument('-i', '--index', help='Column index of source data (one-based)', required=True, type=int)
args = parser.parse_args()

decbins = ['1-9', '10-99', '100-999', '1k+', '10k+', '100k+', '1M+', '10M+', '100M+']

def main():
    bincounts = [0] * len(decbins)
    for line in sys.stdin:
        if line[0] == '#':
            continue
        line = line[:-1]
        parts = line.split()
        val = float(parts[args.index - 1])
        bini = int(math.log10(val))
        if bini < 0 or bini >= len(decbins):
            print('Error: unexpected value {} (acceptable 1.0 -- {}-1'.format(val, 10**len(decbins)-1), file=sys.stderr)
            return
        bincounts[bini] += 1

    lastnonzero = 0
    for i, v in enumerate(bincounts):
        if v > 0:
            lastnonzero = i

    print('# BinName    Count   BinMin   BinMax')
    for i in range(lastnonzero + 1):
        print('{:<9s} {:8d} {:8d} {:8d}'.format(
            decbins[i], bincounts[i], 10**i, 10**(i+1)-1))

if __name__ == '__main__':
    main()
