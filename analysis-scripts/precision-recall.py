#!/usr/bin/python3
# Derive precision and recall from output of compare-safety-runner

import argparse
import io
import json
import os
import string

parser = argparse.ArgumentParser(
    description='Derive precision and recall from output of compare-safety-runner.py')
parser.add_argument('-i', '--indir', help='Directory containing the JSON output files', required=True)
parser.add_argument('-r', '--refdir', help='Directory containing the JSON reference folding files', required=True)
parser.add_argument('-o', '--outdir', help='Directory for output files', required=True)
args = parser.parse_args()

def main():
    ff = os.listdir(args.indir)
    refnames = os.listdir(args.refdir)
    namemap = make_name_map(ff, refnames)
    basedir, rnatype = os.path.split(args.indir)
    if rnatype == '':
        rnatype = os.path.basename(basedir)

    perfile = []
    tmfe = {'tp': 0, 'tn': 0, 'fp': 0, 'fn': 0, 'safe': 0, 'decisions': 0, 'folds': 0}
    tmfesingle = {'tp': 0, 'tn': 0, 'fp': 0, 'fn': 0, 'safe': 0, 'decisions': 0, 'folds': 0}
    tmaxpairs = {'tp': 0, 'tn': 0, 'fp': 0, 'fn': 0, 'safe': 0, 'decisions': 0, 'folds': 0}
    tmaxpairssingle = {'tp': 0, 'tn': 0, 'fp': 0, 'fn': 0, 'safe': 0, 'decisions': 0, 'folds': 0}
    resources = []
    counts = []
    
    for fbase in ff:
        fname = os.path.join(args.indir, fbase)
        refname = os.path.join(args.refdir, namemap[fbase])
        
        with io.open(fname, 'r', encoding='utf-8') as infile, io.open(refname, 'r', encoding='utf-8') as reffile:
            d = json.load(infile)
            reffold = json.load(reffile)
            if not reffold['ReferenceFolding']:
                print('Error: {} lacks a reference folding'.format(fbase))
                continue
            mfe = hitcounts(d.get('RNAsubopt'), reffold)
            mfesingle = singlehitcounts(d.get('RNAsuboptSingle'), reffold)
            maxpairs = hitcounts(d.get('SafeComplete'), reffold)
            maxpairssingle = pairarrayhitcounts(d.get('SingleMaxPairs'), reffold)
            perfile.append({
                'name': reffold['Name'],
                'MFE': mfe, 'MFESingle': mfesingle,
                'MaxPairs': maxpairs, 'MaxPairsSingle': maxpairssingle})
            addto(tmfe, mfe)
            addto(tmfesingle, mfesingle)
            addto(tmaxpairs, maxpairs)
            addto(tmaxpairssingle, maxpairssingle)
            res = getres(d)
            res['Name'] = reffold['Name']
            resources.append(res)
            c = getcounts(d)
            c['Name'] = reffold['Name']
            counts.append(c)

    print('Total: safety as predictor of decision correctness')
    print('                     Precision   Recall SafeDeci NumFoldings TruePositive TrueNegative FalsePositive FalseNegative')
    print('Minimum free energy: {:9.4f} {:8.4f} {:8.4f} {:11d} {:12d} {:12d} {:13d} {:13d}'.format(
        partofsum(tmfe['tp'], tmfe['fp']),
        partofsum(tmfe['tp'], tmfe['fn']),
        divorzero(tmfe['safe'], tmfe['decisions']),
        tmfe['folds'], tmfe['tp'], tmfe['tn'], tmfe['fp'], tmfe['fn']))
    print('Single MFE:          {:9.4f} {:8.4f} {:8.4f} {:11d} {:12d} {:12d} {:13d} {:13d}'.format(
        partofsum(tmfesingle['tp'], tmfesingle['fp']),
        partofsum(tmfesingle['tp'], tmfesingle['fn']),
        divorzero(tmfesingle['safe'], tmfesingle['decisions']),
        tmfesingle['folds'], tmfesingle['tp'], tmfesingle['tn'], tmfesingle['fp'], tmfesingle['fn']))
    print('Maximum pairs:       {:9.4f} {:8.4f} {:8.4f} {:11d} {:12d} {:12d} {:13d} {:13d}'.format(
        partofsum(tmaxpairs['tp'], tmaxpairs['fp']),
        partofsum(tmaxpairs['tp'], tmaxpairs['fn']),
        divorzero(tmaxpairs['safe'], tmaxpairs['decisions']),
        tmaxpairs['folds'], tmaxpairs['tp'], tmaxpairs['tn'], tmaxpairs['fp'], tmaxpairs['fn']))
    print('Max pairs single:    {:9.4f} {:8.4f} {:8.4f} {:11d} {:12d} {:12d} {:13d} {:13d}'.format(
        partofsum(tmaxpairssingle['tp'], tmaxpairssingle['fp']),
        partofsum(tmaxpairssingle['tp'], tmaxpairssingle['fn']),
        divorzero(tmaxpairssingle['safe'], tmaxpairssingle['decisions']),
        tmaxpairssingle['folds'], tmaxpairssingle['tp'], tmaxpairssingle['tn'], tmaxpairssingle['fp'], tmaxpairssingle['fn']))

    print()
    print('Resources taken (average)')
    print('                      CPU-User  CPU-Sys  Mem(kB)')
    print('ViennaRNA RNAsubopt:  {:8.2f} {:8.2f} {:8.0f}'.format(
        avg(resources, 'RNAsuboptUser'),
        avg(resources, 'RNAsuboptSys'),
        avg(resources, 'RNAsuboptMem')))
    print('Trivial safety:       {:8.2f} {:8.2f} {:8.0f}'.format(
        avg(resources, 'TrivialUser'),
        avg(resources, 'TrivialSys'),
        avg(resources, 'TrivialMem')))
    print('Safe & Complete:      {:8.2f} {:8.2f} {:8.0f}'.format(
        avg(resources, 'SCUser'),
        avg(resources, 'SCSys'),
        avg(resources, 'SCMem')))
    print('Single Maximum Pairs: {:8.2f} {:8.2f} {:8.0f}'.format(
        avg(resources, 'SingleUser'),
        avg(resources, 'SingleSys'),
        avg(resources, 'SingleMem')))

    print()
    print('Average counts:   Bases MaxPairs RNAsuboptFoldings SafeCompleteFoldings')
    print('                {:7.1f} {:8.1f} {:17.1f} {:20.1f}'.format(
        avg(counts, 'Bases'),
        avg(counts, 'MaxPairs'),
        avg(counts, 'MFEFolds'),
        avg(counts, 'SCFolds')))

    write_precision_recall(perfile, os.path.join(args.outdir, rnatype+'-precision-recall-mfe.dat'), 'MFE')
    write_precision_recall(perfile, os.path.join(args.outdir, rnatype+'-precision-recall-mfe-single.dat'), 'MFESingle')
    write_precision_recall(perfile, os.path.join(args.outdir, rnatype+'-precision-recall-maxpairs.dat'), 'MaxPairs')
    write_precision_recall(perfile, os.path.join(args.outdir, rnatype+'-precision-recall-maxpairs-single.dat'), 'MaxPairsSingle')
    write_resources(resources, os.path.join(args.outdir, rnatype+'-resources.dat'))
    write_counts(counts, os.path.join(args.outdir, rnatype+'-counts.dat'))

class NAFormatter(string.Formatter):
    def format_field(self, value, format_spec):
        if value == None:
            l = format_spec[:-1].split('.')[0]
            return format('NA', '>' + str(l) + 's')
        try:
            r = format(value, format_spec)
            return r
        except TypeError as e:
            print('Failed to write {} with format {}: {}'.format(value, format_spec, e))
        return 'Error!'

def avg(data, name):
    a = [x.get(name) for x in data]
    nn = [x for x in a if x != None]
    if len(nn) == 0:
        print('No values to average for', name)
        return 0
    return sum(nn) / len(nn)

def write_precision_recall(data, fname, t):
    fmt = NAFormatter()
    with io.open(fname, 'w', encoding='utf-8') as of:
        tee(of, '# Name       Precision   Recall SafeDeci NumFoldings TruePositive TrueNegative FalsePositive FalseNegative')
        for d in data:
              tee(of, fmt.format('{:<12s} {:9.4f} {:8.4f} {:8.4f} {:11d} {:12d} {:12d} {:13d} {:13d}',
                  d['name'],
                  partofsum(d[t]['tp'], d[t]['fp']),
                  partofsum(d[t]['tp'], d[t]['fn']),
                  divorzero(d[t]['safe'], d[t]['decisions']),
                  d[t]['folds'], d[t]['tp'], d[t]['tn'], d[t]['fp'], d[t]['fn']))

def tee(f, s):
    #print(s)
    f.write(s + '\n')

def divorzero(a, b):
    if a == None or b == None:
        if a != None or b != None:
            print("Nullity conflict in division: {} / {}".format(a, b))
        return None
    if b == 0:
        if a != 0:
            print("Error: tried to divide {} / {}".format(a, b)) 
        return 0
    return a / b

def partofsum(a, b):
    if a == None or b == None:
        if a != None or b != None:
            print("Nullity conflict in partofsum: a={}, b={}".format(a, b))
        return None
    if a == 0:
        return 0
    return a / (a + b)

def make_name_map(a, b):
    abase = {}
    for n in a:
        abase[n.split('.')[0]] = n
    bbase = {}
    for n in b:
        bbase[n.split('.')[0]] = n

    ret = {}
    for n in abase:
        ret[abase[n]] = bbase[n]
    return ret

def hitcounts(fold, ref):
    if fold == None:
        return {'tp': None, 'tn': None, 'fp': None, 'fn': None, 'safe': None, 'decisions': None, 'folds': None}
    c = {'tp': 0, 'tn': 0, 'fp': 0, 'fn': 0, 'safe': 0, 'decisions': 0, 'folds': fold['NumFolds']}
    if fold['Bases'] != len(ref['ReferenceFolding']) or fold['Bases'] != len(fold['Free']):
        print('Error: {} has {} bases, reference folding of {} bases and free array of {} bases'.format(
            ref['Name'], fold['Bases'], len(ref['ReferenceFolding']), len(fold['Free'])))
        return {'tp': None, 'tn': None, 'fp': None, 'fn': None, 'safe': None, 'decisions': None, 'folds': None}
    for i in range(fold['Bases']):
        for j in range(i+1, fold['Bases']):
            n = fold['Pairs'][i][j]
            c['decisions'] += n
            if n == fold['NumFolds']:
                c['safe'] += n
                if ref['ReferenceFolding'][i] == j:
                    c['tp'] += n
                else:
                    c['fp'] += n
            else:
                if ref['ReferenceFolding'][i] == j:
                    c['fn'] += n
                else:
                    c['tn'] += n
        n = fold['Free'][i]
        c['decisions'] += n
        if n == fold['NumFolds']:
            c['safe'] += n
            if ref['ReferenceFolding'][i] < 0:
                c['tp'] += n
            else:
                c['fp'] += n
        else:
            if ref['ReferenceFolding'][i] < 0:
                c['fn'] += n
            else:
                c['tn'] += n
    return c

def singlehitcounts(fold, ref):
    if fold == None:
        return {'tp': None, 'tn': None, 'fp': None, 'fn': None, 'safe': None, 'decisions': None, 'folds': None}
    if fold['Bases'] != len(ref['ReferenceFolding']) or fold['Bases'] != len(fold['Free']):
        print('Error: {} has {} bases, reference folding of {} bases and free array of {} bases'.format(
            ref['Name'], fold['Bases'], len(ref['ReferenceFolding']), len(fold['Free'])))
        return {'tp': None, 'tn': None, 'fp': None, 'fn': None, 'safe': None, 'decisions': None, 'folds': None}
    numcorr = 0
    numbad = 0
    preddecis = 0
    for i in range(fold['Bases']):
        for j in range(i+1, fold['Bases']):
            n = fold['Pairs'][i][j]
            preddecis += n
            if ref['ReferenceFolding'][i] == j:
                numcorr += n
            else:
                numbad += n
        n = fold['Free'][i]
        preddecis += n
        if ref['ReferenceFolding'][i] < 0:
            numcorr += n
        else:
            numbad += n
    refdecis = sum([(1 if j < 0 or i < j else 0) for i, j in enumerate(ref['ReferenceFolding'])])
    return {'tp': numcorr, 'tn': 0, 'fp': numbad, 'fn': refdecis - numcorr,
            'safe': preddecis, 'decisions': preddecis, 'folds': fold['NumFolds']}

def pairarrayhitcounts(fold, ref):
    if fold == None:
        return {'tp': None, 'tn': None, 'fp': None, 'fn': None, 'safe': None, 'decisions': None, 'folds': None}
    if fold['Bases'] != len(ref['ReferenceFolding']) or fold['Bases'] != len(fold['Sample']):
        print('Error: {} has {} bases, reference folding of {} bases and sample folding of {} bases'.format(
            ref['Name'], fold['Bases'], len(ref['ReferenceFolding']), len(fold['Sample'])))
        return {'tp': None, 'tn': None, 'fp': None, 'fn': None, 'safe': None, 'decisions': None, 'folds': None}
    numcorr = 0
    numbad = 0
    for i, j in enumerate(fold['Sample']):
        if j < 0:
            if ref['ReferenceFolding'][i] < 0:
                numcorr += 1
            else:
                numbad += 1
        elif i < j:
            if ref['ReferenceFolding'][i] == j:
                numcorr += 1
            else:
                numbad += 1
            
    preddecis = sum([(1 if j < 0 or i < j else 0) for i, j in enumerate(fold['Sample'])])
    refdecis = sum([(1 if j < 0 or i < j else 0) for i, j in enumerate(ref['ReferenceFolding'])])
    return {'tp': numcorr, 'tn': 0, 'fp': numbad, 'fn': refdecis - numcorr,
            'safe': preddecis, 'decisions': preddecis, 'folds': 1}

def addto(total, part):
    for n in total:
        if part[n] != None:
            total[n] += part[n]

def dg(d, *args):
    for k in args:
        if d != None:
            d = d.get(k)
    return d

def getres(d):
    return {
        'RNAsuboptUser': dg(d, 'RNAsubopt', 'Resources', 'RNAsuboptUser'),
        'RNAsuboptSys': dg(d, 'RNAsubopt', 'Resources', 'RNAsuboptSys'),
        'RNAsuboptMem': dg(d, 'RNAsubopt', 'Resources', 'RNAsuboptRSS'),
        'RNAsuboptNum': dg(d, 'RNAsubopt', 'NumFolds'),
        'TrivialUser': dg(d, 'RNAsubopt', 'Resources', 'TrivialSafetyUser'),
        'TrivialSys': dg(d, 'RNAsubopt', 'Resources', 'TrivialSafetySys'),
        'TrivialMem': dg(d, 'RNAsubopt', 'Resources', 'TrivialSafetyRSS'),
        'SCUser': dg(d, 'SafeComplete', 'Resources', 'User'),
        'SCSys': dg(d, 'SafeComplete', 'Resources', 'Sys'),
        'SCMem': dg(d, 'SafeComplete', 'Resources', 'RSS'),
        'SCNum': dg(d, 'SafeComplete', 'NumFolds'),
        'SingleUser': dg(d, 'SingleMaxPairs', 'Resources', 'User'),
        'SingleSys': dg(d, 'SingleMaxPairs', 'Resources', 'Sys'),
        'SingleMem': dg(d, 'SingleMaxPairs', 'Resources', 'RSS'),
    }

def write_resources(res, fname):
    fmt = NAFormatter()
    with io.open(fname, 'w', encoding='utf-8') as of:
        tee(of, '# Name       RNAsuboptUser RNAsuboptSys RNAsuboptMem RNAsuboptNum TrivialUser TrivialSys TrivialMem  SCUser   SCSys   SCMem    SCNum SingleUser SingleSys SingleMem')
        for r in res:
            tee(of, fmt.format('{:<12s} {:13.2f} {:12.2f} {:12d} {:12d} {:11.2f} {:10.2f} {:10d} {:8.2f} {:8.2f} {:8d} {:8d} {:10.2f} {:9.2f} {:9d}',
                r['Name'],
                r['RNAsuboptUser'], r['RNAsuboptSys'], r['RNAsuboptMem'], r['RNAsuboptNum'],
                r['TrivialUser'], r['TrivialSys'], r['TrivialMem'],
                r['SCUser'], r['SCSys'], r['SCMem'], r['SCNum'],
                r['SingleUser'], r['SingleSys'], r['SingleMem']))

def getcounts(d):
    return {
        'Bases': dg(d, 'SafeComplete', 'Bases'),
        'MaxPairs': dg(d, 'SafeComplete', 'NumPairs'),
        'SCFolds': dg(d, 'SafeComplete', 'NumFolds'),
        'MFEFolds': dg(d, 'RNAsubopt', 'NumFolds'),
    }

def write_counts(counts, fname):
    fmt = NAFormatter()
    with io.open(fname, 'w', encoding='utf-8') as of:
        tee(of, '# Name       Bases MaxPairs RNAsuboptFoldings SafeCompleteFoldings')
        for c in counts:
            tee(of, fmt.format('{:<12s} {:5d} {:8d} {:17d} {:20d}',
                               c['Name'], c['Bases'], c['MaxPairs'], c['MFEFolds'], c['SCFolds']))
        
if __name__ == '__main__':
    main()
