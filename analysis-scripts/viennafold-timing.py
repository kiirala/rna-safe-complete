#!/usr/bin/python3
# Run ViennaRNA RNAsubopt and time the run

import argparse
import io
import os
import psutil
import subprocess
import sys
import time

parser = argparse.ArgumentParser(
    description='Run RNAsubopt from ViennaRNA package for number of structures and time results')
parser.add_argument('-r', '--rnasubopt', default='RNAsubopt', help='Location of RNAsubopt program')
parser.add_argument('-i', '--indir', help='Input directory where FASTA files are', required=True)
parser.add_argument('-o', '--outdir', help='Output directory', required=True)
parser.add_argument('-e', '--deltaenergy', help='deltaEnergy parameter to pass to RNAsubopt', type=int, default=1)
args = parser.parse_args()

def countlines(fname):
    count = 0
    with io.open(fname, 'r', encoding='utf-8') as f:
        while f.readline():
            count += 1
    return count

def main():
    fastanames = os.listdir(args.indir)
    print('Analyzing {} sequences'.format(len(fastanames)))

    result = []
    for name in fastanames:
        #print("Analyzing {}... ".format(name), end='')
        ifname = os.path.join(args.indir, name)
        ofname = os.path.join(args.outdir, name.replace('.dp', '.dat'))

        p_args = ['nice', args.rnasubopt, '-e', str(args.deltaenergy), '-i', ifname]
        print(' '.join(p_args), end=' ')
        sys.stdout.flush()

        with io.open(ofname, 'w') as of:
            start = time.time()
            p = psutil.Popen(p_args, stdout=of, stderr=sys.stderr)
            pid, status, res = os.wait4(p.pid, 0)
            elapsed = time.time() - start
            p.wait()
        print('Status: {}, user: {:5.1f}, sys: {:5.1f}, maxrss: {:6d}'.format(
            status, res.ru_utime, res.ru_stime, res.ru_maxrss))
        result.append((name, res.ru_utime, res.ru_stime, res.ru_utime + res.ru_stime, elapsed,
                       res.ru_maxrss, countlines(ofname) - 2))
        #break

    with io.open(os.path.join(args.outdir, 'stats.dat'), 'w', encoding='utf-8') as f:
        f.write('# Name          UserTime SystTime  CPUTime WallTime   MaxRSS Solutions\n')
        for r in result:
            f.write('{:<15s} {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8d} {:9d}\n'.format(*r))
        
if __name__ == '__main__':
    main()
