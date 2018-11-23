#!/usr/bin/python3
# Run ViennaRNA RNAsubopt and time the run

import argparse
import io
import json
import multiprocessing
import os
import psutil
import resource
import signal
import subprocess
import sys
import time

parser = argparse.ArgumentParser(
    description='Runs RNAsubopt from ViennaRNA and keltainen.duckdns.org/rnafolding/comparesafety and stores their timing, memory usage and statistics of their output')
parser.add_argument('-r', '--rnasubopt', default='RNAsubopt', help='Location of ViennaRNA RNAsubopt program')
parser.add_argument('-s', '--safecomplete', default='comparesafety', help='Location of keltainen.duckdns.org/rnafolding/comparesafety program')
parser.add_argument('-t', '--trivialsafety', default='trivialsafety', help='Location of keltainen.duckdns.org/rnafolding/trivialsafety program')
parser.add_argument('-i', '--indir', help='Input directory where FASTA files are', required=True)
parser.add_argument('-o', '--outdir', help='Output directory', required=True)
parser.add_argument('-e', '--deltaenergy', help='deltaEnergy parameter to pass to RNAsubopt', type=int, default=1)
parser.add_argument('-n', '--shard', help='Shard number and total number of shards, e.g. 2:8 for shard 2 out of 8')
parser.add_argument('-c', '--clean', help='Start from clean state, without reading any existing results', type=bool)
parser.add_argument('-w', '--workers', help='Number of simultaneous worker processes to start', type=int, default=1)
parser.add_argument('-u', '--timeout', help='Timeout for RNAsubopt in hours. In case of timeout, the sequence is skipped', type=float, default=0)
args = parser.parse_args()

def countlines(fname):
    count = 0
    with io.open(fname, 'r', encoding='utf-8') as f:
        while f.readline():
            count += 1
    return count

def run_rnasubopt(ifname, deltaenergy, number=None):
    errs = []
    rna_args = ['nice', args.rnasubopt, '-e', str(deltaenergy), '-i', ifname]
    safety_args = ['nice', args.trivialsafety]
    if number != None:
        safety_args.append('-num')
        safety_args.append(str(number))

    #print(' '.join(rna_args), '|', ' '.join(safety_args))

    rna = psutil.Popen(rna_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    safety = psutil.Popen(safety_args, stdin=rna.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    rna.stdout.close()

    if args.timeout > 0:
        resource.prlimit(rna.pid, resource.RLIMIT_CPU, (int(args.timeout*60*60), int(args.timeout*60*60)))

    folddata = None
    try:
        folddata = json.loads(safety.stdout.read().decode('utf-8'))
    except json.decoder.JSONDecodeError as e:
        errs.append('{}: Failed to decode trivial safety output:'.format(ifname, e))

    rnaerr = rna.stderr.read()
    if rnaerr != None and len(rnaerr) > 0:
        errs.append('RNAsubopt for {} returned errors:\n{}'.format(ifname, rnaerr))
    safetyerr = safety.stderr.read()
    if safetyerr != None and len(safetyerr) > 0:
        errs.append('Trivial safety for {} returned errors:\n{}'.format(ifname, safetyerr))

    rpid, rstatus, rres = os.wait4(rna.pid, 0)
    #print('RNA:    Status: {}, user (s): {:5.1f}, sys (s): {:5.1f}, maxrss (kB): {:6d}'.format(
    #    rstatus, rres.ru_utime, rres.ru_stime, rres.ru_maxrss))
    spid, sstatus, sres = os.wait4(safety.pid, 0)
    #print('Safety: Status: {}, user (s): {:5.1f}, sys (s): {:5.1f}, maxrss (kB): {:6d}'.format(
    #    sstatus, sres.ru_utime, sres.ru_stime, sres.ru_maxrss))

    if os.WIFSIGNALED(rstatus) and not (number != None and os.WTERMSIG(rstatus) == signal.SIGPIPE):
        errs.append('{}: RNAsubopt was terminated with signal {}'.format(ifname, os.WTERMSIG(rstatus)))
        return (None, errs)
    if os.WIFSIGNALED(sstatus):
        errs.append('{}: Trivialsafety was terminated with signal {}'.format(ifname, os.WTERMSIG(sstatus)))
        return (None, errs)
    if folddata == None:
        errs.append('{}: folddata was None for unknown reason'.format(ifname))
        return (None, errs)

    folddata['Command'] = ' '.join(rna_args), '|', ' '.join(safety_args)
    folddata['Resources'] = {
        'RNAsuboptUser': rres.ru_utime,
        'RNAsuboptSys': rres.ru_stime,
        'RNAsuboptRSS': rres.ru_maxrss,
        'TrivialSafetyUser': sres.ru_utime,
        'TrivialSafetySys': sres.ru_stime,
        'TrivialSafetyRSS': sres.ru_maxrss,
    }
    return (folddata, errs)

def run_safecomplete(ifname):
    errs = []
    sc_args = ['nice', args.safecomplete, '-json', '-in', ifname]
    #print(' '.join(sc_args))

    sc = psutil.Popen(sc_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    folddata = None
    try:
        folddata = json.loads(sc.stdout.read().decode('utf-8'))
    except json.decoder.JSONDecodeError as e:
        errs.append('{}: Failed to decode trivial safety output:'.format(ifname, e))

    err = sc.stderr.read()
    if err != None and len(err) > 0:
        errs.append('Safe&Complete for {} returned errors:\n{}'.format(ifname, err))

    pid, status, res = os.wait4(sc.pid, 0)
    #print('SaCmpl: Status: {}, user (s): {:5.1f}, sys (s): {:5.1f}, maxrss (kB): {:6d}'.format(
    #    status, res.ru_utime, res.ru_stime, res.ru_maxrss))

    if folddata == None:
        return (None, errs)

    folddata['Command'] = ' '.join(sc_args)
    folddata['Resources'] = {
        'User': res.ru_utime,
        'Sys': res.ru_stime,
        'RSS': res.ru_maxrss,
    }
    return (folddata, errs)
    
def run_singlemaxpairs(ifname):
    errs = []
    sc_args = ['nice', args.safecomplete, '-json', '-single', '-in', ifname]
    #print(' '.join(sc_args))

    sc = psutil.Popen(sc_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    folddata = None
    try:
        folddata = json.loads(sc.stdout.read().decode('utf-8'))
    except json.decoder.JSONDecodeError as e:
        errs.append('{}: Failed to decode trivial safety output:'.format(ifname, e))

    err = sc.stderr.read()
    if err != None and len(err) > 0:
        errs.append('Single-Max-Pairs for {} returned errors:\n{}'.format(ifname, err))

    pid, status, res = os.wait4(sc.pid, 0)
    #print('SMP:    Status: {}, user (s): {:5.1f}, sys (s): {:5.1f}, maxrss (kB): {:6d}'.format(
    #    status, res.ru_utime, res.ru_stime, res.ru_maxrss))

    if folddata == None:
        return (None, errs)

    folddata['Command'] = ' '.join(sc_args)
    folddata['Resources'] = {
        'User': res.ru_utime,
        'Sys': res.ru_stime,
        'RSS': res.ru_maxrss,
    }
    return (folddata, errs)

def run_single(name):
    errs = []
    ifname = os.path.join(args.indir, name)
    ofname = os.path.join(args.outdir, name.replace('.dp', '.json').replace('.fasta', '.json'))

    vf, svf, sc, smp = (None, None, None, None)
    if not args.clean:
        of = None
        try:
            of = io.open(ofname, 'r', encoding='utf-8')
        except FileNotFoundError:
            pass #errs.append('File {} doesn\'t exist, will run analysis'.format(ofname))
        except Exception as e:
            errs.append('Failed to open file {}: {}'.format(ofname, e))
        if of:
            try:
                d = json.load(of)
                vf = d.get('RNAsubopt')
                svf = d.get('RNAsuboptSingle')
                sc = d.get('SafeComplete')
                smp = d.get('SingleMaxPairs')
            except Exception as e:
                errs.append('Failed to read existing data from {}: {}'.format(ofname, e))

    has = [1 if x != None else 0 for x in (vf, sc, smp)]
    if sum(has) > 0 and sum(has) < 3:
        errs.append('Partial data exists: Viennafold: {}, Safe&Complete: {}, Single-Max-Pairs: {}'.format(*has))
    if not vf:
        vf, vfe = run_rnasubopt(ifname, args.deltaenergy)
        errs += vfe
    if not svf:
        svf, vfe = run_rnasubopt(ifname, 0, number=1)
        errs += vfe
    if not sc:
        sc, sce = run_safecomplete(ifname)
        errs += sce        
    if not smp:
        smp, smpe = run_singlemaxpairs(ifname)
        errs += smpe
        
    if vf != None and vf.get('Bases') != sc.get('Bases'):
        errs.append('ViennaRNA found {} bases, Safe&Complete {} bases!'.format(vf.get('Bases'), sc.get('Bases')))
    if vf != None and vf.get('Name') != sc.get('Name'):
        errs.append('ViennaRNA found name {}, Safe&Complete name {}!'.format(vf.get('Name'), sc.get('Name')))

    with io.open(ofname, 'w', encoding='utf-8') as f:
        json.dump({'RNAsubopt': vf, 'RNAsuboptSingle': svf, 'SafeComplete': sc, 'SingleMaxPairs': smp}, f)

    return name, vf, sc, smp, errs

def main():
    fastanames = os.listdir(args.indir)
    if not args.shard:
        print('Analyzing {} sequences'.format(len(fastanames)), file=sys.stderr)
    else:
        ss = args.shard.split(':')
        if len(ss) != 2 or not ss[0].isdigit() or not ss[1].isdigit():
            print('Shard number and count should look like 2:8 (for shard two out of eight)', file=sys.stderr)
            return
        shardi, shardn = [int(s) for s in ss]
        if shardi < 0 or shardi > shardn:
            print('Shard number should be between 0 and shard count', file=sys.stderr)
        if shardi == shardn:
            shardi = 0
        fastalen = len(fastanames)
        fastanames = fastanames[shardi::shardn]
        print('Analyzing {} out of {} sequences (shard {} of {})'.format(len(fastanames), fastalen, shardi, shardn), file=sys.stderr)
    sys.stderr.flush()

    print('# Name        VF_Secs VF_RSSkB  TS_Secs TS_RSSkB  SC_Secs SC_RSSkB SMP_Secs SMP_RSSk     VF_folds     SC_folds')
    sys.stdout.flush()
    
    result = []
    with multiprocessing.Pool(processes=args.workers) as pool:
        for name, vf, sc, smp, errs in pool.imap_unordered(run_single, fastanames):
            if len(errs) > 0:
                print('{} produced {} errors:\n{}'.format(name, len(errs), "\n".join(errs)), file=sys.stderr)
            sys.stdout.flush()
            sys.stderr.flush()
            if vf == None or sc == None or smp == None:
                continue
            print('{:<12s} {:8.1f} {:8d} {:8.1f} {:8d} {:8.1f} {:8d} {:8.1f} {:8d} {:12d} {:12d}'.format(
                vf['Name'],
                vf['Resources']['RNAsuboptUser'] + vf['Resources']['RNAsuboptSys'],
                vf['Resources']['RNAsuboptRSS'],
                vf['Resources']['TrivialSafetyUser'] + vf['Resources']['TrivialSafetySys'],
                vf['Resources']['TrivialSafetyRSS'],
                sc['Resources']['User'] + sc['Resources']['Sys'],
                sc['Resources']['RSS'],
                smp['Resources']['User'] + smp['Resources']['Sys'],
                smp['Resources']['RSS'],
                vf['NumFolds'],
                sc['NumFolds'],
            ))
            sys.stdout.flush()
            sys.stderr.flush()

        #break
        
if __name__ == '__main__':
    main()
