#!/usr/bin/python3
import argparse
import io
import json
import os
import sys

parser = argparse.ArgumentParser(
    description='Analyze locations of safe bases within tRNA foldings')
parser.add_argument('-o', '--outdir', help='Write results to files in named directory', required=True)
parser.add_argument('-s', '--safecomplete', help='Directory containing JSON output files of folding.go', required=True)
parser.add_argument('-m', '--mfe', help='Directory containing MFE foldings from ViennaFold RNAsubopt program', required=True)
args = parser.parse_args()

def main():
    scstats = []
    mfestats = []
    for fname in os.listdir(args.safecomplete):
        scfname = os.path.join(args.safecomplete, fname)
        mfefname = os.path.join(args.mfe, fname.replace('.json', '.fasta'))
        with open(scfname, 'r', encoding='utf-8') as scf, open(mfefname, 'r', encoding='utf-8') as mfef:
            sc = json.load(scf)
            sca = analyze([x['Pairing'] for x in sc['AllFoldings']],
                          sc['ReferenceFolding']['Pairing'], sc['ReferencePosition'])
            addto(scstats, sca)
            mfe = readRnaSubopt(mfef)
            mfea = analyze(mfe.folds, sc['ReferenceFolding']['Pairing'], sc['ReferencePosition'])
            addto(mfestats, mfea)

    with io.open(os.path.join(args.outdir, 'safe-locations.dat'), 'w', encoding='utf-8') as of:
        tee(of, '# Base NumSeqs MPNumSafe MPSafeCorrect MPUnsafeCorrect MFENumSafe MFESafeCorrect MFEUnsafeCorrect')
        for i, s in enumerate(zip(scstats, mfestats)):
            tee(of, '{:<6d} {:7d} {:9d} {:13.2f} {:15.2f} {:10d} {:14.62} {:16.2f}'.format(i, *s[0], *s[1][1:]))
        
def tee(f, text):
    print(text)
    f.write(text + "\n")

def samepairing(a, b, pos):
    if a[pos] < 0 and b[pos] < 0:
        return True
    if a[pos] == b[pos]:
        return True
    return False

def analyze(foldings, reference, positions):
    res = []
    for i in range(max(positions) + 1):
        pos = -1
        found = False
        issafe = False
        safecount = 0
        safecorr = 0
        unsafecorr = 0
        try:
            pos = positions.index(i)
        except ValueError:
            pass
        if pos >= 0:
            found = True
            safecount = sum([1 if samepairing(x, foldings[0], pos) else 0 for x in foldings])
            issafe = safecount == len(foldings)
            corr = sum([1 if samepairing(x, reference, pos) else 0 for x in foldings])
            if issafe:
                safecorr = corr
            else:
                unsafecorr = corr
        res.append([1 if found else 0,
                    1 if issafe else 0,
                    safecorr / len(foldings),
                    unsafecorr / len(foldings)])
    return res

def addto(total, entry):
    while len(total) < len(entry):
        total.append([0] * len(entry[0]))
    for i, v in enumerate(entry):
        for j in range(len(v)):
            total[i][j] += v[j]

class SubOpt:
    def __init__(self, comment, e, seq, folds):
        self.comment = comment
        self.e = e
        self.seq = seq
        self.folds = folds

def dotBracketToIndices(s):
    leftpos = []
    ans = [-1] * len(s)
    for i, c in enumerate(s):
        if c == '(':
            leftpos.append(i)
        elif c == ')':
            j = leftpos.pop()
            ans[i] = j
            ans[j] = i
        elif c == '.':
            pass
        else:
            print('dotBracketToIndices: unexpected character "{}"'.format(c))
    if len(leftpos) != 0:
        print('dotBracketToIndices: values left in stack at end')
    return ans

def readRnaSubopt(f):
    ll = [x.strip() for x in f.readlines()]
    if ll[0][0] != '>':
        print('readRnaSubopt: expectation failed: comment doesn\'t start with > character')
    comment = ll[0][1:]
    seq, _, e = ll[1].split()
    e = float(e)
    folds = [dotBracketToIndices(l.split()[0]) for l in ll[2:]]
    for fold in folds:
        if len(seq) != len(fold):
            print('readRnaSubopt: expectation failed: sequence and folding have different lengths')
    return SubOpt(comment, e, seq, folds)


if __name__ == '__main__':
    main()
