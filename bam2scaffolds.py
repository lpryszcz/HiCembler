#!/usr/bin/env python
desc="""Report scaffolds by joining contigs based on contact matrix from BAM file. 

TBD:
- allow for fitting contig into gap within large contig?
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Bratislava, 27/10/2016
"""

# Force matplotlib to not use any Xwindows backend.
import matplotlib; matplotlib.use('Agg')

import ete3, glob, gzip, os, resource, subprocess, sys
import scipy.cluster.hierarchy as sch
import fastcluster
import numpy as np
from collections import Counter
from datetime import datetime
from multiprocessing import Pool
from FastaIndex import FastaIndex
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from bam2clusters import logger, normalize, array2tree, bam2clusters, bam2arrays, \
     estimate_distance_parameters, contact_func, distance_func

def normalize_window_size(d, bin_chr, bin_position, windowSize=0):
    """Return symmetric and normalised matrix by window size"""
    logger("Window-size normalisation...")
    # make symmetric & normalise
    d += d.T
    d -= np.diag(d.diagonal()/2)
    #return d, bin_chr, bin_position
    # normalize by windows size
    sizes = np.diff(bin_position, axis=1)#[:, 0]
    if not windowSize:
        c = Counter(sizes.reshape(len(sizes)))
        windowSize, occurencies = c.most_common(1)[0]; print windowSize, occurencies
    d *= 1. * windowSize **2 / (sizes*sizes.T)
    return d, bin_chr, bin_position    

def normalize_rows(a):
    """Normalise rows so the sums among rows are identical."""
    rows, cols = a.shape
    maxv = a.sum(axis=0).max()
    for i in xrange(rows):
        # only if any signal
        if a[i].max():
            a[i] *= 1.*maxv/a[i].sum()
    return a
            
def _func_reduce(A, keys, func, allkeys=None):
    """Reduces along first dimension by aggregating rows with the same keys.
    new row order will be sorted by keys,  i.e. given by: np.unique(keys)
    """
    unique_keys = np.unique(keys)
    if allkeys == None:
        allkeys = unique_keys
    newshape = (len(allkeys), ) + A.shape[1:]
    newA = np.zeros(newshape, dtype = A.dtype)
    for i, k in enumerate(allkeys):
        indices  =  (keys == k)       
        newA[i] = func(A[indices], axis = 0)
    return newA
        
def _average_reduce(A, keys):
    return _func_reduce(A, keys, func=np.mean)
    
def _average_reduce_2d(A, keys):
    return _average_reduce(_average_reduce(A, keys).T, keys).T
        
def fasta2windows(fasta, windowSize, verbose, skipShorter=1, minSize=2000,
                  contigs=[], filterwindows=1):
    """Generate windows over chromosomes"""
    # init fasta index
    faidx = FastaIndex(fasta)
    if verbose:
        logger("Parsing FastA file...")
    # filter windows so they are smaller than largest chr and withing reasonalbe range toward genome size
    if filterwindows:    
        maxchrlen = max(faidx.id2stats[c][0] for c in faidx)
        windowSize = filter(lambda x: 1000*x<maxchrlen and 1000*x<0.01*faidx.genomeSize and 1000*x>0.000001*faidx.genomeSize, windowSize)
        if verbose:
            logger(" selected %s windows [kb]: %s"%(len(windowSize), str(windowSize)))
    windowSize = [w*1000 for w in windowSize]
    # generate windows
    windows, chr2window = [[] for w in windowSize], [{} for w in windowSize]
    genomeSize = 0
    base2chr = {}
    skipped = []
    for i, c in enumerate(faidx, 1): 
        if contigs and c not in contigs:
            continue
        if i%1e5 == 1:
            sys.stderr.write(' %s   \r'%i)
        # get windows
        size = faidx.id2stats[c][0]
        # skip short contigs    
        if skipShorter and size < min(windowSize) or size<minSize:
            skipped.append(size)
            continue
        for ii in range(len(windowSize)):
            # skip contig if shorter than given window
            if skipShorter and size < windowSize[ii]:
                #print size, "skipped %s"%windowSize[i]
                continue
            # get starting window
            chr2window[ii][c] = len(windows[ii])
            for start in range(0, size, windowSize[ii]):
                windows[ii].append((c, start, start+windowSize[ii]))
            # update last entry end
            if not skipShorter:
                windows[ii][-1] = (c, start, size)
        # get chromosome tick in the middle    
        base2chr[genomeSize+size/2] = c
        # update genomeSize
        genomeSize += size
    # print [len(w) for w in windows]
    if verbose:
        logger(' %s bases in %s contigs divided in %s-%s windows. '%(faidx.genomeSize, len(faidx), len(windows[0]), len(windows[-1])))
        if skipped:
            logger('  %s bases in %s contigs skipped.'%(sum(skipped), len(skipped)))
    return windowSize, windows, chr2window, base2chr, faidx.genomeSize

def get_reversed(scaffold):
    """Return reversed scaffold, updating orientation"""
    return [(name, not orientation) for name, orientation in reversed(scaffold)]

def get_indices(scaffold, contig2indices):
    """Return list with indices representing given scaffold"""
    indices = []
    for name, reverse in scaffold:
        if name not in contig2indices:
            continue
        if reverse:
            indices += reversed(contig2indices[name])
        else:
            indices += contig2indices[name]
    return indices
    
def join_scaffolds(scaffold1, scaffold2, d, contig2indices, minWindows=3):
    """Join two adjacent scaffolds"""
    indices1 = get_indices(scaffold1, contig2indices)
    indices2 = get_indices(scaffold2, contig2indices)
    # skip contigs with less windows than minWindows
    if len(indices1) < len(indices2):
        scaffold1, indices1, scaffold2, indices2 = scaffold2, indices2, scaffold1, indices1
    # skip scaffolding if len of contig smaller than necessary
    if len(indices2) < minWindows:
        if len(indices1) < minWindows:
            return []
        return scaffold1
    # get subset of array for scaffold1 and scaffold2
    _d = d[:, indices1+indices2][indices1+indices2, :]
    # get orientation: 0: s-s; 1: s-e; 2: e-s; 3: e-e
    n1, n2 = len(indices1), len(indices2)
    # compare sums of subsets of contact matrix
    i = n2/2
    d1, d2 = _d[:n1/2, n1:], _d[n1-n1/2:n1, n1:]
    orientation = np.argmax(map(np.sum, (d1[:,:i], d1[:,-i:], d2[:,:i], d2[:,-i:])))
    # s - s
    if   orientation == 0:
        scaffold = get_reversed(scaffold2) + scaffold1 
    # s - e
    elif orientation == 1:
        scaffold = scaffold2 + scaffold1 
    # e - s
    elif orientation == 2:
        scaffold = scaffold1 + scaffold2 
    # e - e
    else:
        scaffold = scaffold1 + get_reversed(scaffold2)
    return scaffold

def get_contig2indices(bin_chr):
    """Return contig2 indices"""
    contig2indices = {c: [] for c in np.unique(bin_chr)}
    for i, c in enumerate(bin_chr):
        contig2indices[c].append(i)
    return contig2indices
    
def tree2scaffold(t, d, bin_chr, bin_position, minWindows):
    """Scaffold contigs based on distance matrix and tree."""    
    # get contig with indices
    contig2indices = get_contig2indices(bin_chr)

    # populate internal nodes with growing scaffolds
    for n in t.traverse('postorder'):
        # add scaffold for each leave
        # each scaffold consists of ordered list of contigs and their orientations (0/False: Fwd or 1/True:Rev)
        if n.is_leaf():
            n.scaffold = [(n.name, 0)]
            continue
        # unload children
        n1, n2 = n.get_children()
        # and combine scaffolds
        n.scaffold = join_scaffolds(n1.scaffold, n2.scaffold, d, contig2indices, minWindows)

    return t.scaffold
        
def scaffold2gaps(d, bin_chr, scaffold, params, fasta, mingap=100, maxgap=1000000, gapfrac=.33):
    """Add gaps to scaffolds
    Here, there may be some smarter way of doing it using ML and all contig contacts at the same time. 
    """
    faidx = FastaIndex(fasta)
    contig2size = {c: stats[0] for c, stats in faidx.id2stats.iteritems()}
    contig2indices = get_contig2indices(bin_chr)
    scaffold_with_gaps = []
    contigbp = gapbp = 0
    for i, (c, o) in enumerate(scaffold[:-1]):
        c2, o2 = scaffold[i+1]
        csize, c2size = contig2size[c], contig2size[c2]
        contacts = d[contig2indices[c][0],contig2indices[c2][0]]
        distance = distance_func(contacts, *params)*1000
        gapSize = distance - contig2size[c]/2 - contig2size[c2]/2
        if gapSize < mingap:
            gapSize = mingap
        elif gapSize > gapfrac * (contig2size[c] + contig2size[c2]):
            gapSize = gapfrac * (contig2size[c] + contig2size[c2])
        if gapSize > maxgap:
            gapSize = maxgap
        # round to int
        gapSize = int(round(gapSize))
        contigbp += csize
        gapbp += gapSize
        scaffold_with_gaps.append((c, o, gapSize))#; print gapSize, csize, c2size, distance
    # add last contig without gap
    c, o = scaffold[-1]
    scaffold_with_gaps.append((c, o, 0))
    #sys.exit()
    return scaffold_with_gaps, contigbp, gapbp
    
def contigs2scaffold(args):
    """Combine contigs into scaffold"""
    i, bam, fasta, windowSize, contigs, mapq, minWindows, params = args
    
    # get windows
    data = fasta2windows(fasta, windowSize, verbose=0, skipShorter=1, contigs=set(contigs), filterwindows=0)
    windowSize, windows, chr2window, base2chr, genomeSize = data
    
    # get array from bam
    arrays = [np.zeros((len(w), len(w)), dtype="float32") for w in windows]
    arrays = bam2arrays(arrays, windowSize, chr2window, bam,  mapq, regions=contigs, verbose=0)
    d = arrays[0]
    
    # generate missing handles
    bin_chr, bin_position = [], []
    for c, s, e in windows[0]:
        bin_chr.append(c)
        bin_position.append((s, e))
    bin_chr = np.array(bin_chr)
    bin_position = np.array(bin_position)
    
    # skip if empty array or report without scaffolding if only one contig
    if not d.shape[0]:
        return []
    if len(np.unique(bin_chr)) == 1:
        return [(bin_chr[0], 0, 0)]
        
    # make symmetric & normalise
    d += d.T - np.diag(d.diagonal())

    # get tree on reduced matrix
    transform = lambda x: distance_func(x+1, *params)
    t = array2tree(transform(_average_reduce_2d(d, bin_chr)), np.unique(bin_chr))
    
    # get scaffold
    scaffold = tree2scaffold(t, d, bin_chr, bin_position, minWindows)
    
    # estimate & add gaps
    #scaffold = [(c, o, 1) for c, o in scaffold]; return scaffold
    scaffold, contigbp, gapbp = scaffold2gaps(_average_reduce_2d(d, bin_chr), np.unique(bin_chr), scaffold, params, fasta)
    
    info = " cluster_%s with %s windows for %s out of %s contigs; %s bp in gaps [%2.2f%s]"
    logger(info%(i, d.shape[0], len(np.unique(bin_chr)), len(contigs), gapbp, 100.*gapbp/contigbp, '%'))

    return scaffold

def get_scaffolds(out, bam, fasta, clusters, mapq, minWindows, threads, verbose, windowSize=10):
    """Return scaffolds"""
    scaffolds = []
    faidx = FastaIndex(fasta)
    contig2size = {c: stats[0] for c, stats in faidx.id2stats.iteritems()}
    
    # estimate distance parameters
    params = estimate_distance_parameters(out, bam, mapq, contig2size, windowSize*1000)
    # 
    if threads > 1:
        p = Pool(threads)
        iterable = [(i, bam, fasta, [windowSize], contigs, mapq, minWindows, params) for i, contigs in enumerate(clusters, 1)]
        for scaffold in p.imap(contigs2scaffold, iterable):
            scaffolds.append(scaffold)
    else:
        for i, contigs in enumerate(clusters, 1): 
            scaffold = contigs2scaffold([i, bam, fasta, [windowSize], contigs, mapq, minWindows, params])
            scaffolds.append(scaffold)
    return scaffolds

def report_scaffolds(outbase, scaffolds, faidx, w=60):
    """Save scaffolds"""
    totsize = 0
    fastafn = outbase+".scaffolds.fa"
    scaffoldfn = outbase+".scaffolds.tab"
    with open(fastafn, "w") as out, open(scaffoldfn, "w") as outscaffolds:
        for i, scaffold in enumerate(scaffolds, 1):
            # skip empty scaffolds
            if not scaffold:
                continue
            seqs = []
            elements = []
            for c, o, gapSize in scaffold:
                seqs.append(faidx.get_sequence(c, reverse=o))
                seqs.append('N'*gapSize)
                if o:
                    elements.append("%s-"%c)
                else:
                    elements.append("%s+"%c)
            # store seq and architecture
            seq = "".join(seqs)
            seq = "\n".join(seq[s:s+w] for s in range(0, len(seq), w))
            out.write(">scaffold%s %s bp in %s contigs\n%s\n"%(i, len(seq), len(elements), seq))
            outscaffolds.write("%s\n"%"\t".join(elements))
            totsize += len(seq)
    logger(" %s in %s scaffolds reported to %s"%(totsize, len(scaffolds), fastafn))
    return fastafn
    
def bam2scaffolds(bam, fasta, outdir, minSize, windowSizes, mapq, threads, dpi, upto, minWindows, verbose):
    """Report scaffolds based on BAM file"""
    faidx = FastaIndex(fasta)

    # check windows size
    maxwindow = max(stats[0] for stats in faidx.id2stats.itervalues())/ 10e3
    windowSizes = filter(lambda x: x <= maxwindow, windowSizes)
    if not windowSizes:
        info = "Select window sizes (-w) <= %s kb (at least 10x smaller than the longest contig)!\n"
        sys.stderr.write(info%int(maxwindow))
        sys.exit(1)
        
    logger(" Selected %s window sizes <= %s kb: %s"%(len(windowSizes), int(maxwindow), str(windowSizes)))

    # calculate clusters for various window size
    clusters = bam2clusters(bam, fasta, outdir, minSize, mapq, threads, dpi, upto, verbose)
    
    for windowSize in windowSizes:
        outbase = os.path.join(outdir,"auto.%sk"%windowSize)
        
        logger("=== %sk windows ==="%windowSize)
        logger("Constructing scaffolds...")
        scaffolds = get_scaffolds(outbase, bam, fasta.name, clusters, mapq, minWindows, threads, verbose, windowSize)

        logger("Reporting %s scaffolds..."%len(scaffolds))
        fastafn = report_scaffolds(outbase, scaffolds, faidx)
    
    logger("Done!")

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--version', action='version', version='1.01e')   
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose")    
    parser.add_argument("-i", "--bam", nargs="+", help="BAM file(s)")
    parser.add_argument("-f", "--fasta", type=file, help="Genome FastA file")
    parser.add_argument("-o", "--outdir", required=1, help="output name")
    parser.add_argument("-w", "--windowSize", default=[5, 10, 2], type=int, nargs="+", 
                        help="window size in kb used for scaffolding [%(default)s]")
    parser.add_argument("-m", "--minSize", default=2000, type=int,
                        help="minimum contig length [%(default)s]")
    parser.add_argument("-q", "--mapq", default=10, type=int,
                        help="mapping quality [%(default)s]")
    parser.add_argument("-u", "--upto", default=0, type=float,
                        help="process up to this number of reads from each library [all]")
    parser.add_argument("-t", "--threads", default=4, type=int,
                        help="number of processes to use [%(default)s]")
    parser.add_argument("-d", "--dpi", default=300, type=int,
                        help="output images dpi [%(default)s]")
    parser.add_argument("--minWindows", default=2, type=int,
                        help="minimum number of windows per contig, has to be >1 [%(default)s]")

    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    # check bam files
    missing = [_bam for _bam in o.bam if not os.path.isfile(_bam)]
    if missing:
        sys.stderr.write("No such file(s): %s\n"%", ".join(missing))
        sys.exit(1)
        
    # create outdir
    if os.path.isdir(o.outdir):
        sys.stderr.write("Output directory exists: %s !\n"%o.outdir)
        #sys.exit(1)
    else:
        os.makedirs(o.outdir)
        
    # process
    bam2scaffolds(o.bam, o.fasta, o.outdir, o.minSize, o.windowSize, o.mapq, o.threads,
                  o.dpi, o.upto, o.minWindows, o.verbose)

if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    #except IOError as e:
    #    sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror))
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
    