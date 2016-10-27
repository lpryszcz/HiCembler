#!/usr/bin/env python
desc="""Report scaffolds by joining contigs based on contact matrix from BAM file. 

TBD:
- add SAM support
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Bratislava, 27/10/2016
"""

import ete3, glob, gzip, os, resource, subprocess, sys
import scipy.cluster.hierarchy as sch
import numpy as np
from collections import Counter
from datetime import datetime
from multiprocessing import Pool
from fastq2array import logger, fasta2windows, get_window, plot
from array2scaffolds import transform, normalize_rows, get_name, get_clusters, array2tree, tree2scaffold, _average_reduce_2d, report_scaffolds
from FastaIndex import FastaIndex

def _get_samtools_proc(bam, mapq=0, regions=[], skipFlag=3980):
    """Return samtools subprocess"""
    # skip second in pair, unmapped, secondary, QC fails and supplementary algs
    args = map(str, ['samtools', 'view', "-q", mapq, "-F", skipFlag, bam])
    # add regions
    if regions:
        args += regions
    # start subprocess
    proc = subprocess.Popen(args, stdout=subprocess.PIPE)
    return proc

def bam2array(arrays, windowSize, chr2window, bam, mapq, upto=1e7, regions=[], verbose=1):
    """Return contact matrix based on BAM"""
    i = 1
    # init empty array
    for _bam in bam:
        # init samtools
        if verbose:
            logger(" %s"%_bam)
        proc = _get_samtools_proc(_bam, mapq, regions)
        # process reads from given library
        for i, line in enumerate(proc.stdout, i):
            if upto and i>upto:
                break
            if verbose and not i%1e5:
                sys.stderr.write(" %s \r"%i)
            # unload data
            ref1, start1, mapq, cigar, ref2, start2, insertSize, seq = line.split('\t')[2:10]
            start1, start2, seqlen = int(start1), int(start2), len(seq)
            # update ref if alignment in the same chrom
            if ref2 == "=":
                ref2 = ref1
            # skip if contig not present in array
            if ref1 not in chr2window[0] or ref2 not in chr2window[0]:
                continue
            # get windows
            for ii in range(len(windowSize)):
                w1 = get_window(ref1, start1, seqlen, windowSize[ii], chr2window[ii])
                w2 = get_window(ref2, start2, seqlen, windowSize[ii], chr2window[ii])
                # matrix is symmetrix, so make sure to add only to one part
                if w2 < w1:
                    w1, w2 = w2, w1
                # update contact array
                arrays[ii][w1][w2] += 1
    if verbose:
        logger(" %s alignments parsed"%i)
    # stop subprocess
    proc.terminate()
    return arrays

def load_clusters(fnames):
    """Load clusters from directory"""
    clusters, windowSize = [], []
    for fn in fnames:
        clusters.append([l[:-1].split('\t') for l in open(fn)])
        w = int(os.path.basename(fn).split('k')[0])*1000
        windowSize.append(w)
    return clusters, windowSize
    
def get_clusters(bam, fasta, outdir, windowSize, mapq, dpi, upto, verbose):
    """Return clusters computed from from windowSizes"""
    clusters = []
    # load clusters
    fnames = glob.glob(os.path.join(outdir, "*k.clusters.tab"))
    if fnames:
        clusters, windowSize = load_clusters(fnames)
        return clusters, windowSize
    
    # get windows
    windowSize, windows, chr2window, base2chr, genomeSize = fasta2windows(fasta, windowSize, verbose, skipShorter=0)

    # get array from bam
    logger("Parsing BAM...")
    arrays = [np.zeros((len(w), len(w)), dtype="float32") for w in windows]
    arrays = bam2array(arrays, windowSize, chr2window, bam,  mapq, upto)

    faidx = FastaIndex(fasta)
    contig2size = {c: stats[0] for c, stats in faidx.id2stats.iteritems()}
    for d, _windowSize, _windows in zip(arrays, windowSize, windows):
        outfn = outdir + "/%sk"%(_windowSize/1000,)
        logger("=== %sk windows ==="%(_windowSize/1000,))
        '''# save windows, array and plot
        logger("Saving & plotting...")
        with gzip.open(outfn+".windows.tab.gz", "w") as out:
            out.write("\n".join("\t".join(map(str, w)) for w in _windows)+"\n")
        with open(outfn+".npz", "w") as out:
            np.savez_compressed(out, d)
        if len(windows)<2e4:
            plot(outfn, a, genomeSize, base2chr, _windowSize, dpi)
        elif verbose:
            sys.stderr.write("[WARNING] Very large matrix (%s x %s). Skipped plotting!\n"%(len(windows), len(windows)))
        #'''
        
        # make symmetric & normalise
        d += d.T
        d -= np.diag(d.diagonal()/2)
        d = normalize_rows(d)

        # generate missing handles
        bin_chr, bin_position = [], []
        for c, s, e in _windows:
            bin_chr.append(c)
            bin_position.append((s, e))
        bin_chr = np.array(bin_chr)
        bin_position = np.array(bin_position)

        logger("Calculating linkage matrix & tree...")
        names = ["%s %7sk"%(get_name(c), s/1000) for c, (s, e) in zip(bin_chr, bin_position)]
        d = transform(d)
        t = array2tree(d, names, outfn)
        del d

        logger("Assigning contigs to clusters/scaffolds...")
        _clusters = get_clusters(outfn, t, contig2size, bin_chr)
        clusters.append(_clusters)
    return clusters, windowSize

def contigs2scaffold(args):
    """Combine contigs into scaffold"""
    i, bam, fasta, windowSize, contigs, mapq, minWindows = args
    # get array from bam
    # get windows
    data = fasta2windows(fasta, windowSize, verbose=0, skipShorter=0, contigs=set(contigs), filterwindows=0)
    windowSize, windows, chr2window, base2chr, genomeSize = data

    #logger("  Parsing BAM...")
    arrays = [np.zeros((len(w), len(w)), dtype="float32") for w in windows]
    arrays = bam2array(arrays, windowSize, chr2window, bam,  mapq, regions=contigs, verbose=0)
    d = arrays[0]
    # make symmetric & normalise
    d += d.T
    d -= np.diag(d.diagonal()/2)
    d = normalize_rows(d)

    logger(" cluster_%s with %s windows in %s contigs"%(i, d.shape[0], len(contigs)))
    # generate missing handles
    bin_chr, bin_position = [], []
    for c, s, e in windows[0]:
        bin_chr.append(c)
        bin_position.append((s, e))
    bin_chr = np.array(bin_chr)
    bin_position = np.array(bin_position)

    # get tree on reduced matrix
    t = array2tree(transform(_average_reduce_2d(d, bin_chr)), np.unique(bin_chr))

    # get scaffold
    scaffold = tree2scaffold(t, d, bin_chr, bin_position, minWindows)
    return scaffold

def get_scaffolds(bam, fasta, clusters, mapq, minWindows, threads, verbose, windowSize=[10]):
    """Return scaffolds"""
    scaffolds = []
    faidx = FastaIndex(fasta)
    contig2size = {c: stats[0] for c, stats in faidx.id2stats.iteritems()}
    if threads>1:
        p = Pool(threads)
        iterable = [(i, bam, fasta, windowSize, contigs, mapq, minWindows) for i, contigs in enumerate(clusters, 1)]
        for scaffold in p.imap(contigs2scaffold, iterable):
            scaffolds.append(scaffold)
    else:
        for i, contigs in enumerate(clusters, 1):
            scaffold = contigs2scaffold([i, bam, fasta, windowSize, contigs, mapq, minWindows])
            scaffolds.append(scaffold)
    return scaffolds
    
def bam2scaffolds(bam, fasta, outdir, windowSize, windowSize2, mapq, threads, dpi, upto, minWindows, verbose):
    """Report scaffolds based on BAM file"""
    faidx = FastaIndex(fasta)

    # calculate clusters for various window size
    clusters, _windowSize = get_clusters(bam, fasta, outdir, windowSize, mapq, dpi, upto, verbose)

    # use preselected window size
    if len(windowSize)==1 and windowSize[0]*1000 in _windowSize:
        selected = _windowSize.index(windowSize[0]*1000)
    else:
        selected = np.argmin(map(len, clusters)) # check also number of contigs & cummulative size
        
    # choose window with least scaffolds
    logger("=== Scaffolding ===")
    _clusters, _windowSize = clusters[selected], _windowSize[selected]
    logger("Selected %s clusters from window %sk"%(len(_clusters), _windowSize/1000))

    logger("Constructing scaffolds...")
    scaffolds = get_scaffolds(bam, fasta.name, _clusters, mapq, minWindows, threads, verbose, [windowSize2])
    
    logger("Reporting %s scaffolds..."%len(scaffolds))
    outbase = outdir+"/%sk.%sk"%(_windowSize/1000, windowSize2)
    fastafn = report_scaffolds(outbase, scaffolds, faidx)
    
    logger("Done!")

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--version', action='version', version='1.0b')   
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose")    
    parser.add_argument("-i", "--bam", nargs="+", help="BAM file(s)")
    parser.add_argument("-f", "--fasta", type=file, help="Genome FastA file")
    parser.add_argument("-o", "--outdir", required=1, help="output name")
    parser.add_argument("-w", "--windowSize", nargs="+", default=[100, 50, 20, 10], type=int,
                        help="window size in kb [%(default)s]")
    parser.add_argument("-z", "--windowSize2", default=10, type=int,
                        help="window size in kb [%(default)s]")
    parser.add_argument("-m", "--mapq", default=10, type=int,
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
        
    # create outdir
    if os.path.isdir(o.outdir):
        sys.stderr.write("Output directory exists: %s !\n"%o.outdir)
        #sys.exit(1)
    else:
        os.makedirs(o.outdir)
        
    # process
    bam2scaffolds(o.bam, o.fasta, o.outdir, o.windowSize, o.windowSize2, o.mapq, o.threads, o.dpi, o.upto, o.minWindows, o.verbose)

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
    