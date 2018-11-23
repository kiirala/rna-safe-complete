#!/usr/bin/python3
import argparse
import io
import json
import numpy
import os
import sys

parser = argparse.ArgumentParser(
    description='Analyze predicted safety or RNA secondary structure vs. biological ground truth')
parser.add_argument('-o', '--outdir', help='Write results to files in named directory', required=True)
parser.add_argument('files', help='JSON files with RNA secondary structure to analyze', nargs='+')
args = parser.parse_args()

def tee(f, text):
    print(text)
    f.write(text + "\n")

def main():
    results = []
    of = io.open(os.path.join(args.outdir, 'stats.dat'), 'w', encoding='utf-8')
    tee(of, '# Name NumFoldings SafeBases BaseCorr SafeCorr UnsafeCorr')
    for fname in args.files:
        with io.open(fname, 'r', encoding='utf-8') as f:
            data = json.load(f)
            r = analyze(data)
            tee(of, '%6s %11d %9.6f %8.6f %8.6f %10.6f' % (
                r['Name'],
                data['Counts']['SafeCompleteFoldings'],
                float(data['Counts']['SafeBases']) / data['Counts']['SequenceBases'],
                numpy.average(r['PBaseCorrect']),
                r['SafeCorrectness'], r['UnsafeCorrectness']))
            results.append(r)
    of.close()
    
    with io.open(os.path.join(args.outdir, 'dump.json'), 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2)

    #summ = summarize(results)
    #with io.open(os.path.join(args.outdir, 'basesafety.dat'), 'w', encoding='utf-8') as f:
    #    pass

def printperc(a):
    print(' '.join(['%3.0f' % (x * 100) for x in a]))

def pairequal(a, b):
    return a == b or (a < 0 and b < 0)

def partof(a, b):
    if a > 0 or b > 0:
        return a / (a + b)
    return 0.0

def analyze(d):
    ref = d['ReferenceFolding']['Pairing']
    pbase = [0.0] * len(ref)
    pfree = []
    nnpfree = []
    ppair = []
    safecorrect = [[0.0, 0.0], [0.0, 0.0]]
    freesafecorrect = []
    pairsafecorrect = []
    fractocorrect = []
    
    for i in range(len(ref)):
        if ref[i] < 0:
            frac = float(d['Safety']['FreeCount'][i]) / d['Counts']['SafeCompleteFoldings']
            pbase[i] = frac
            pfree.append(frac)
            if d['Sequence'][i] != 'N':
                nnpfree.append(frac)
                freesafecorrect.append([d['Safety']['SafeBase'][i], frac])
        else:
            posa = min(i, ref[i])
            posb = max(i, ref[i])
            frac = float(d['Safety']['PairCount'][posa][posb]) / d['Counts']['SafeCompleteFoldings']
            pbase[i] = frac
            if i < ref[i]:
                ppair.append(frac)
                pairsafecorrect.append([d['Safety']['SafeBase'][i], frac])

        if d['Sequence'][i] == 'N':
            pass
        elif d['Safety']['SafeBase'][i]:
            safecorrect[0][0] += pbase[i]
            safecorrect[0][1] += 1 - pbase[i]
        else:
            safecorrect[1][0] += pbase[i]
            safecorrect[1][1] += 1 - pbase[i]

    for i in range(len(ref)):
        ffrac = float(d['Safety']['FreeCount'][i]) / d['Counts']['SafeCompleteFoldings']
        if ffrac > 0:
            fractocorrect.append((ffrac, ref[i] < 0, d['Sequence'][i] == 'N'))
        for j in range(i+1, len(ref)):
            pfrac = float(d['Safety']['PairCount'][i][j]) / d['Counts']['SafeCompleteFoldings']
            if pfrac > 0:
                fractocorrect.append((pfrac, ref[i] == j, False))

    #print(numpy.average(pbase)*100, numpy.average(ppair)*100, numpy.average(pfree)*100, numpy.average(nnpfree)*100)
    #print(safecorrect,
    #      partof(safecorrect[0][0], safecorrect[0][1]),
    #      partof(safecorrect[1][0], safecorrect[1][1]))

    return {'Name': d['Name'],
            'PBaseCorrect': pbase,
            'PFreeCorrect': pfree,
            'PNonNFreeCorrect': nnpfree,
            'PPairCorrect': ppair,
            'SafeCorrectCount': safecorrect,
            'SafeCorrectness': partof(safecorrect[0][0], safecorrect[0][1]),
            'UnsafeCorrectness': partof(safecorrect[1][0], safecorrect[1][1]),
            'FreeSafeCorrect': freesafecorrect,
            'PairSafeCorrect': pairsafecorrect,
            'FracToCorrect': fractocorrect}
    
if __name__ == '__main__':
    main()

