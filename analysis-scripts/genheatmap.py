#!/usr/bin/python3

import argparse
import io
import json
import sys

parser = argparse.ArgumentParser(
    description='Generate a GNUplot heatmap data from pairing arrays in JSON format')
parser.add_argument('-i', '--input', help='Input file', required=True)
parser.add_argument('-o', '--output', help='Output file', required=True)
parser.add_argument('-p', '--path', help='Dot-separated path of the pairing array within the input file', required=True)
parser.add_argument('-r', '--reference', help='Sequence and reference folding file in JSON format', required=True)
args = parser.parse_args()

bases = ['-', 'A', 'C', 'G', 'U', 'N']

def main():
    with io.open(args.input, 'r', encoding='utf-8') as infile, io.open(args.reference, 'r', encoding='utf-8') as reffile:
        d = json.load(infile)
        ref = json.load(reffile)
        arr = d
        for n in args.path.split('.'):
            if n:
                arr = arr.get(n)
        if not arr:
            print('Path {} not found in file {}', args.path, args.input, file=sys.stdout)
            return
        for i, n in enumerate(arr.get('Free')):
            arr.get('Pairs')[i][i] = n
        if not ref.get('Bases'):
            print('Array \'Bases\' not found in file {}', args.input, file=sys.stdout)
            
        with io.open(args.output, 'w', encoding='utf-8') as outfile:
            writearray(arr.get('Pairs'), ref.get('Bases'), outfile)

def writearray(arr, ref, out):
    refc = ['{{\\\\tiny{{}}{}}}'.format(bases[x]) for x in ref]
    out.write('X {}\n'.format(' '.join(refc)))
    for base, line in zip(refc, arr):
        out.write('{} {}\n'.format(base, ' '.join([str(x) if x else '0' for x in line])))

if __name__ == '__main__':
    main()
