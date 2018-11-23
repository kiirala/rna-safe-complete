#!/usr/bin/python3
# Compute safety from the output of ViennaFold RNAsubopt output.

import argparse
import io
import json
import os

parser = argparse.ArgumentParser(
    description='Compute safety from the output of ViennaFold RNAsubopt output')
parser.add_argument('-s', '--seq', help='Directory containing sequence files', required=True)
parser.add_argument('-f', '--fold', help='Directory containing directories of folding files', required=True)
parser.add_argument('-o', '--outdir', help='Write results to files in named directory', required=True)
args = parser.parse_args()

def main():
    fdirs = os.listdir(args.fold)
    seqfiles = os.listdir(args.seq)
    print('Analyzing {} sequences'.format(len(seqfiles)))

    res = []
    for basename in seqfiles:
        foldfnames = [os.path.join(args.fold, fd, basename) for fd in fdirs]
        seqfname = os.path.join(args.seq, basename)
        with io.open(seqfname, 'r', encoding='utf-8') as fseq:
            seq = readFasta(fseq)
            for foldfname in foldfnames:
                with io.open(foldfname, 'r', encoding='utf-8') as ffold:
                    fold = readRnaSubopt(ffold)
                    res.append(analyze(seq, fold))

    with io.open(os.path.join(args.outdir, 'viennafold-safety-dump.json'), 'w', encoding='utf-8') as of:
        json.dump(res, of, indent=2, default=lambda o: o.__dict__)

    # Number of correct predictions found at certain predicted percentage
    with io.open(os.path.join(args.outdir, 'viennafold-found-correct.dat'), 'w', encoding='utf-8') as of:
        of.write('# P CorrAtE1 CorrBlv1 CorrAtE2 CorrBlv2 CorrAtE3 CorrBlv3 CorrAtE4 CorrBlv4 CorrAtE5 CorrBlv5\n')
        for p in range(0, 101):
            at = [frac_correct_between(res, a, p / 100.0, (p+1) / 100.0) for a in [1, 2, 3, 4, 5]]
            blv = [frac_correct_between(res, a, 0, (p+1) / 100.0) for a in [1, 2, 3, 4, 5]]
            of.write('{:<3d} {:8.6f} {:8.6f} {:8.6f} {:8.6f} {:8.6f} {:8.6f} {:8.6f} {:8.6f} {:8.6f} {:8.6f}\n'.format(
                p, at[0], blv[0], at[1], blv[1], at[2], blv[2], at[3], blv[3], at[4], blv[4]))

    # Number of correct predictions out of all predictions at certain predicted percentage
    with io.open(os.path.join(args.outdir, 'viennafold-frac-correct.dat'), 'w', encoding='utf-8') as of:
        of.write('# P CorrAtE1 CorrBlv1 CorrAtE2 CorrBlv2 CorrAtE3 CorrBlv3 CorrAtE4 CorrBlv4 CorrAtE5 CorrBlv5\n')
        for p in range(0, 101):
            at = [frac_correct_predictions_of_all(res, a, p / 100.0, (p+1) / 100.0) for a in [1, 2, 3, 4, 5]]
            blv = [frac_correct_predictions_of_all(res, a, 0, (p+1) / 100.0) for a in [1, 2, 3, 4, 5]]
            of.write('{:<3d} {:8.6f} {:8.6f} {:8.6f} {:8.6f} {:8.6f} {:8.6f} {:8.6f} {:8.6f} {:8.6f} {:8.6f}\n'.format(
                p, at[0], blv[0], at[1], blv[1], at[2], blv[2], at[3], blv[3], at[4], blv[4]))

    with io.open(os.path.join(args.outdir, 'frac-safe.dat'), 'w', encoding='utf-8') as of:
        of.write('# E SafeBases SafePredictions\n')
        for e in [1, 2, 3, 4, 5]:
            safebases, safepreds = frac_safe(res, e)
            of.write('{:<3d} {:9.6f} {:15.6f}\n'.format(e, safebases, safepreds))

class Fasta:
    def __init__(self, comment, seq, fold):
        self.comment = comment
        self.seq = seq
        self.fold = fold

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

def readFasta(f):
    ll = [x.strip() for x in f.readlines()]
    if len(ll) % 2 != 1:
        print('readFasta: expectation failed: file has even number of lines')
    if ll[0][0] != '>':
        print('readFasta: expectation failed: comment doesn\'t start with > character')
    comment = ll[0][1:]
    seq = ''.join(ll[1:int(len(ll)/2+1)])
    fold = dotBracketToIndices(''.join(ll[int(len(ll)/2+1):]))
    if len(seq) != len(fold):
        print('readFasta: expectation failed: sequence and its folding are of different length')
    return Fasta(comment, seq, fold)

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

class Analysis:
    def __init__(self, preds, e, bases, decisions):
        self.preds = preds
        self.e = e
        self.bases = bases
        self.decisions = decisions

class Prediction:
    def __init__(self, p, correct, numbases):
        if p < 0 or p > 1:
            print("Created prediction with p={:.2f}".format(p))
        self.p = p
        self.correct = correct
        self.numbases = numbases

def analyze(seq, subopt):
    decisions = sum([1 if i < 0 or i < seq.fold[i] else 0 for i in seq.fold])

    free = [0] * len(seq.seq)
    pair = [[0] * len(seq.seq) for _ in range(len(seq.seq))]
    for fold in subopt.folds:
        for i, j in enumerate(fold):
            if j < 0:
                free[i] += 1
            elif i < j:
                pair[i][j] += 1
    preds = []
    for i in range(len(seq.seq)):
        if free[i] > 0:
            preds.append(Prediction(
                float(free[i]) / len(subopt.folds),
                seq.fold[i] < 0, 1))
        for j in range(i+1, len(seq.seq)):
            if pair[i][j] > 0:
                preds.append(Prediction(
                    float(pair[i][j]) / len(subopt.folds),
                    seq.fold[i] == j, 2))

    return Analysis(preds, subopt.e, len(seq.fold), decisions)

def divorzero(a, b):
    if b == 0:
        return 0
    return float(a) / b

def frac_correct_between(res, e, pmin, pmax):
    ss = sum([r.preds for r in res if r.e == e], [])
    corr = [1 if s.correct else 0 for s in ss if s.p >= pmin and s.p < pmax]
    return float(sum(corr)) / sum([r.decisions for r in res if r.e == e])

def frac_correct_predictions_of_all(res, e, pmin, pmax):
    ss = sum([r.preds for r in res if r.e == e], [])
    corr = [1 if s.correct else 0 for s in ss if s.p >= pmin and s.p < pmax]
    return divorzero(sum(corr), len(corr))

def frac_safe(res, e):
    safebases = 0
    safedecis = 0
    totalbases = 0
    totaldecis = 0
    for r in res:
        if r.e != e:
            continue
        safebases += sum([p.numbases for p in r.preds if p.p == 1.0])
        safedecis += sum([1 for p in r.preds if p.p == 1.0])
        totalbases += r.bases
        totaldecis += r.decisions
    return (safebases / totalbases, safedecis / totaldecis)

if __name__ == '__main__':
    main()
