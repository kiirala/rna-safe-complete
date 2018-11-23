#!/usr/bin/python3
import argparse
import io
import json
import numpy
import os
import random

parser = argparse.ArgumentParser(
    description='Summarize the output from bio-vs-safety.py')
parser.add_argument('dir', help='Directory containing the dump.json output file')
args = parser.parse_args()

def main():
    with io.open(os.path.join(args.dir, 'dump.json'), 'r', encoding='utf-8') as f:
        data = json.load(f)

    predicted_fraction_when_safe(data)
    print()
    print('Fraction of correct predictions having probability between')
    for p in range(0, 101, 5):
        vals = fraction_correct_with_predicted_probability(data, p / 100.0, (p + 5) / 100.0)
        print('{:.2f} {:.2f} {:.4f} {:.4f} {:.4f} {:.4f}'.format(
            p / 100.0, (p + 5) / 100.0, *vals))
    print()
    print('Fraction of correct predictions having probability over')
    for p in range(0, 101, 5):
        vals = fraction_correct_with_predicted_probability(data, p / 100.0, 1.1)
        print('{:.2f} {:.4f} {:.4f} {:.4f} {:.4f}'.format(
            p / 100.0, *vals))

    with io.open(os.path.join(args.dir, 'fraccorrect.dat'), 'w', encoding='utf-8') as f:
        f.write('# P   CorrAt   StdDev NNCorrAt NNStdDev  CorrSum NNCorrSum Shuffled NNShuffled\n')
        for p in range(0, 101):
            at = fraction_correct_with_predicted_probability(data, p / 100.0, (p + 1) / 100.0)
            cumsum = fraction_correct_found_upto_probability(data, p / 100.0)
            scumsum = fraction_correct_found_upto_probability_shuffled(data, p / 100.0)
            f.write('{:<3d} {:8.6f} {:8.6f} {:8.6f} {:8.6f} {:8.6f} {:9.6f} {:8.6f} {:10.6f}\n'.format(
                p, *at, *cumsum, *scumsum))

def predicted_fraction_when_safe(dd):
    fsc = []
    psc = []
    for d in dd:
        fsc += d['FreeSafeCorrect']
        psc += d['PairSafeCorrect']
    cfsc = [x[1] for x in fsc if x[0]]
    ifsc = [x[1] for x in fsc if not x[0]]
    cpsc = [x[1] for x in psc if x[0]]
    ipsc = [x[1] for x in psc if not x[0]]
    print('Fraction of predictions that match biological folding:')
    print('Free base, when predicted safe:     avg={:.4f} stddev={:.4f}'.format(
        numpy.average(cfsc), numpy.std(cfsc)))
    print('Free base, when predicted non-safe: avg={:.4f} stddev={:.4f}'.format(
        numpy.average(ifsc), numpy.std(ifsc)))
    print('Pair, when predicted safe:          avg={:.4f} stddev={:.4f}'.format(
        numpy.average(cpsc), numpy.std(cpsc)))
    print('Pair, when predicted non-safe:      avg={:.4f} stddev={:.4f}'.format(
        numpy.average(ipsc), numpy.std(ipsc)))

def fraction_correct_with_predicted_probability(dd, pmin, pmax):
    ftc = []
    for d in dd:
        ftc += d['FracToCorrect']
    abv = [1 if x[1] else 0 for x in ftc if x[0] >= pmin and x[0] < pmax]
    nnabv = [1 if x[1] else 0 for x in ftc if (x[0] >= pmin and x[0] < pmax) and not x[2]]
    return (numpy.average(abv), numpy.std(abv), numpy.average(nnabv), numpy.std(nnabv))

def fraction_correct_found_upto_probability(dd, pmax):
    ftc = []
    items = 0
    nnitems = 0
    for d in dd:
        ftc += d['FracToCorrect']
        items += len(d['PFreeCorrect']) + len(d['PPairCorrect'])
        nnitems += (len(d['PFreeCorrect']) + len(d['PPairCorrect']) -
                    sum([1 if x[1] and x[2] else 0 for x in d['FracToCorrect']]))
    abv = [1 if x[1] else 0 for x in ftc if x[0] <= pmax]
    nnabv = [1 if x[1] else 0 for x in ftc if x[0] <= pmax and not x[2]]
    return (float(sum(abv)) / items, float(sum(nnabv)) / items)

def fraction_correct_found_upto_probability_shuffled(dd, pmax):
    ftc = []
    items = 0
    nnitems = 0
    for d in dd:
        ftc += d['FracToCorrect']
        items += len(d['PFreeCorrect']) + len(d['PPairCorrect'])
        nnitems += (len(d['PFreeCorrect']) + len(d['PPairCorrect']) -
                    sum([1 if x[1] and x[2] else 0 for x in d['FracToCorrect']]))
    fracshuffled = [x[0] for x in ftc]
    random.shuffle(fracshuffled)
    shuff = [(x, y[1], y[2]) for x, y in zip(fracshuffled, ftc)]
    abv = [1 if x[1] else 0 for x in shuff if x[0] <= pmax]
    nnabv = [1 if x[1] else 0 for x in shuff if x[0] <= pmax and not x[2]]
    return (float(sum(abv)) / items, float(sum(nnabv)) / items)

if __name__ == '__main__':
    main()
