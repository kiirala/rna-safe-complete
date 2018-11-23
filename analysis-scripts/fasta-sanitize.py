#!/usr/bin/python3
# Convert STRAND dot-bracket files to sanitized FASTA format.
# The ViennaRNA package can't read STRAND files directly:
#   * the sequence must be on a single line
#   * only one comment line is allowed

import argparse
import io
import os
import sys

parser = argparse.ArgumentParser(
    description='Sanitize files in FASTA + dot-bracket format')
parser.add_argument('-i', help='Input file', required=True)
parser.add_argument('-o', help='Output file')
args = parser.parse_args()

def main():
    out = sys.stdout
    if args.o:
        out = io.open(args.o, 'w', encoding='utf-8')
    with io.open(args.i, 'r', encoding='utf-8') as inp:
        comment, seq, fold = readFasta(inp, args.i)
        out.write('>{}\n'.format(comment))
        out.write(seq + '\n')
        out.write(fold + '\n')

def readFasta(f, name):
    ll = [x.strip() for x in f.readlines()]
    if ll[0][0] != '>' and ll[0][0] != '#':
        print('readFasta {}: expectation failed: comment doesn\'t start with > or # character'.format(name))

    cml = [x for x in ll if len(x) > 0 and (x[0] == '> ' or x[0] == '#')]
    ncl = [x for x in ll if len(x) > 0 and not (x[0] == '> ' or x[0] == '#')]

    if len(ncl) % 2 != 0:
        print('readFasta {}: expectation failed: file has odd number of non-comment lines'.format(name))

    comment = [x[1:].strip() for x in cml]
    if comment[0].startswith("File "):
        comment[0] = comment[0][5:].strip()
    seq = ''.join(ncl[:len(ncl)//2])
    seq = seq.replace('.', 'N')
    fold = ''.join(ncl[len(ncl)//2:])
    if len(seq) != len(fold):
        print('readFasta {}: expectation failed: sequence and its folding are of different length'.format(name))
    return (' - '.join(comment), seq, fold)

if __name__ == '__main__':
    main()
