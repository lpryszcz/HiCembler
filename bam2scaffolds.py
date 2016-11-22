#!/usr/bin/env python
desc="""Report scaffolds by joining contigs based on contact matrix from BAM file. 

TBD:
- allow for fitting contig into gap within large contig?
- distance estimation?
- remove fastq2array and array2scaffolds imports
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Bratislava, 27/10/2016
"""

import ete3, glob, gzip, os, resource, subprocess, sys
import scipy.cluster.hierarchy as sch
import fastcluster
import numpy as np
from collections import Counter
from datetime import datetime
from multiprocessing import Pool
from fastq2array import logger, fasta2windows, get_window, plot
from array2scaffolds import transform, normalize_rows, get_name, get_clusters, array2tree, _average_reduce_2d, report_scaffolds, getNewick
from FastaIndex import FastaIndex

def normalize_diagional(d, bin_chr, bin_position):
    """Return symmetric and diagonal normalised matrix"""
    #logger("Diagonal normalisation...")    
    # skip rows/columns with zero diagonal
    indices = d.diagonal()!=0
    d = d[indices, :]
    d = d[:, indices]
    bin_chr = bin_chr[indices]
    bin_position = bin_position[indices, :]
    # make symmetric & normalise
    d += d.T
    d -= np.diag(d.diagonal()/2)
    # diag normalisation
    d *= np.mean(d.diagonal()) / d.diagonal()
    return d, bin_chr, bin_position

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

def load_clusters(fnames, _windowSize=[]):
    """Load clusters from directory"""
    clusters, windowSize = [], []
    for fn in fnames:
        w = int(os.path.basename(fn).split('k')[0])*1000
        if _windowSize and w/1000 in _windowSize:
            clusters.append([l[:-1].split('\t') for l in open(fn)])
            windowSize.append(w)
    return clusters, windowSize

def get_longest(t, maxdist=6, k=2.0):
    """Return node having longest branch
    THIS CAN BE FASTER DEFINITELY!
    """
    #n = sorted(t.traverse(), key=lambda n: 2*n.dist-t.get_distance(n), reverse=1)[0]
    n = t
    bestdist = k*n.dist-n.get_distance(t)
    for _n in t.traverse():
        if _n.get_distance(t, topology_only=1) > maxdist:
            break
        if k*_n.dist-_n.get_distance(t) > bestdist:
            n = _n
            bestdist = k*_n.dist-_n.get_distance(t)
    return n, bestdist
            
def get_subtrees(d, bin_chr, bin_position, method="ward", nchrom=1000, distfrac=0.4):
    """Return contings clustered into scaffolds
    fastcluster is slightly faster than scipy.cluster.hierarchy and solve the issue: http://stackoverflow.com/a/40738977/632242
    """
    maxtdist = 0
    i = 0
    subtrees = []
    names = ["%s %s"%(c, s) for c, (s, e) in zip(bin_chr, bin_position)] #names = get_names(bin_chr, bin_position)
    Z = fastcluster.linkage(d[np.triu_indices(d.shape[0], 1)], method=method) 
    # t = distance_matrix2tree(Z, names); print len(set(t.get_leaf_names())), len(names)
    tree = sch.to_tree(Z, False)
    t = ete3.Tree(getNewick(tree, "", tree.dist, names))
    for i in range(1, nchrom):
        tname, tdist = t.get_farthest_leaf()
        if maxtdist < tdist:
            maxtdist = tdist
        # get longest branch
        n, bestdist = get_longest(t)
        # break if small subtree
        if tdist / maxtdist < 1.1 * bestdist / tdist or tdist < maxtdist*distfrac: 
            break
        # store cluster
        subtrees.append(n.get_leaf_names())        
        # prune tree
        # removing child is much faster than pruning and faster than recomputing matrix
        ancestors = n.get_ancestors()
        p = ancestors[0]
        p.remove_child(n)
        n2 = p.get_children()[0]
        if len(ancestors) < 2: 
            p.remove_child(n2)
            t = n2
            t.dist = 0
        else:
            p2 = ancestors[1]
            p2.remove_child(p)
            p2.add_child(n2, dist=n2.dist+p.dist)
    if i:    
        subtrees.append(t.get_leaf_names())
    return subtrees
    
def get_clusters(outbase, d, contig2size, bin_chr, bin_position, method="ward", nchrom=1000, distfrac=0.75):
    """Return clusters for given distance matrix"""
    maxtdist = 0
    d = transform(d)
    subtrees = get_subtrees(d, bin_chr, bin_position, method)

    logger(" Assigning contigs to %s clusters..."%len(subtrees))
    total = correct = 0
    contig2cluster = {get_name(c): Counter() for c in np.unique(bin_chr)}
    for i, subtree in enumerate(subtrees, 1):
        c = Counter(map(get_name, subtree))
        total += len(subtree)
        correct += c.most_common(1)[0][1]
        # poplate contig2cluster
        for k, v in c.iteritems():
            if not k: continue
            contig2cluster[get_name(k)][i] += v
    logger("  %s / %s [%.2f%s]"%(correct, total, 100.*correct/total, '%'))

    logger(" Weak assignments...")
    clusters = [[] for i in range(len(subtrees)+1)]
    withoutCluster, weakCluster = [], []
    for c, counter in contig2cluster.iteritems():
        if not counter:
            withoutCluster.append(c)
            continue
        # get major cluster
        clusteri, count = counter.most_common(1)[0]#; print clusteri, len(clusters), count
        mfrac = 1. * count / sum(counter.itervalues())
        clusters[clusteri].append(c)
        if mfrac < .66:
            weakCluster.append(c)
    logger("  %s bp in %s contigs without assignment."%(sum(contig2size[c] for c in withoutCluster), len(withoutCluster)))
    logger("  %s bp in %s contigs having weak assignment."%(sum(contig2size[c] for c in weakCluster), len(weakCluster)))
      
    outfile = outbase+".clusters.tab"
    clusters = filter(lambda x: x, clusters)
    totsize = 0
    logger("Reporting %s clusters to %s ..."%(len(clusters), outfile))
    with open(outfile, "w") as out:
        for i, cluster in enumerate(clusters, 1):
            #print " cluster_%s %s windows; %s"%(i, len(cluster), Counter(get_chromosome(cluster).most_common(3)))
            clSize = sum(contig2size[c] for c in cluster)
            totsize += clSize
            out.write("\t".join(cluster)+"\n")
    logger("  %3s bp in %s clusters generated from %s contigs."%(totsize, len(clusters), len(contig2cluster)))
    return clusters        
    
def bam2clusters(bam, fasta, outdir, windowSize, mapq, dpi, upto, verbose):
    """Return clusters computed from from windowSizes"""
    clusters = []
    # load clusters
    fnames = glob.glob(os.path.join(outdir, "*k.clusters.tab")) #%",".join(map(str, windowSize))))
    if fnames:
        clusters, windowSize = load_clusters(fnames, windowSize)
        return clusters, windowSize
    
    # get windows
    windowSize, windows, chr2window, base2chr, genomeSize = fasta2windows(fasta, windowSize, verbose, skipShorter=0)
    
    #'''# get array from bam
    logger("Parsing BAM...")
    arrays = [np.zeros((len(w), len(w)), dtype="float32") for w in windows]
    arrays = bam2array(arrays, windowSize, chr2window, bam,  mapq, upto)
    ''' # code for regenerating clusters
    arrays = []
    for _w in windowSize:
        npy = np.load(os.path.join(outdir, "%sk.npz"%(_w/1000,)))
        arrays.append(npy[npy.files[0]])
    print [a.shape for a in arrays]#'''
    
    faidx = FastaIndex(fasta)
    contig2size = {c: stats[0] for c, stats in faidx.id2stats.iteritems()}
    for d, _windowSize, _windows in zip(arrays, windowSize, windows):
        outfn = outdir + "/%sk"%(_windowSize/1000,)
        logger("=== %sk windows ==="%(_windowSize/1000,))
        
        # save windows, array and plot
        logger("Saving array...")
        with gzip.open(outfn+".windows.tab.gz", "w") as out:
            out.write("\n".join("\t".join(map(str, w)) for w in _windows)+"\n")
        with open(outfn+".npz", "w") as out:
            np.savez_compressed(out, d)

        # generate missing handles
        bin_chr, bin_position = [], []
        for c, s, e in _windows:
            bin_chr.append(c)
            bin_position.append((s, e))
        bin_chr = np.array(bin_chr)
        bin_position = np.array(bin_position)

        # make symmetric & normalise
        d, bin_chr, bin_position = normalize_diagional(d, bin_chr, bin_position)

        logger("Assigning contigs to clusters/scaffolds...")
        _clusters = get_clusters(outfn, d, contig2size, bin_chr, bin_position)
        clusters.append(_clusters)
    return clusters, windowSize

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
    if len(indices2) < minWindows:
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

    # estimate distances
    return t.scaffold

def contigs2scaffold(args):
    """Combine contigs into scaffold"""
    i, bam, fasta, windowSize, contigs, mapq, minWindows = args
    # get array from bam
    # get windows
    data = fasta2windows(fasta, windowSize, verbose=0, skipShorter=1, contigs=set(contigs), filterwindows=0)
    windowSize, windows, chr2window, base2chr, genomeSize = data

    #logger("  Parsing BAM...")
    arrays = [np.zeros((len(w), len(w)), dtype="float32") for w in windows]
    arrays = bam2array(arrays, windowSize, chr2window, bam,  mapq, regions=contigs, verbose=0)
    d = arrays[0]
    
    # generate missing handles
    bin_chr, bin_position = [], []
    for c, s, e in windows[0]:
        bin_chr.append(c)
        bin_position.append((s, e))
    bin_chr = np.array(bin_chr)
    bin_position = np.array(bin_position)
    
    # make symmetric & normalise
    d += d.T
    d -= np.diag(d.diagonal()/2)
    d = normalize_rows(d)
    logger(" cluster_%s with %s windows in %s contigs"%(i, d.shape[0], len(contigs)))
    
    # get tree on reduced matrix
    t = array2tree(transform(_average_reduce_2d(d, bin_chr)), np.unique(bin_chr))
    
    # get scaffold
    scaffold = tree2scaffold(t, d, bin_chr, bin_position, minWindows)
    return scaffold

def get_scaffolds(bam, fasta, clusters, mapq, minWindows, threads, verbose, windowSize=10):
    """Return scaffolds"""
    scaffolds = []
    faidx = FastaIndex(fasta)
    contig2size = {c: stats[0] for c, stats in faidx.id2stats.iteritems()}
    if threads > 1:
        p = Pool(threads)
        iterable = [(i, bam, fasta, [windowSize], contigs, mapq, minWindows) for i, contigs in enumerate(clusters, 1)]
        for scaffold in p.imap(contigs2scaffold, iterable):
            scaffolds.append(scaffold)
    else:
        for i, contigs in enumerate(clusters, 1):
            scaffold = contigs2scaffold([i, bam, fasta, [windowSize], contigs, mapq, minWindows])
            scaffolds.append(scaffold)
    return scaffolds
    
def bam2scaffolds(bam, fasta, outdir, windowSize, windowSize2, mapq, threads, dpi, upto, minWindows, verbose):
    """Report scaffolds based on BAM file"""
    faidx = FastaIndex(fasta)

    # calculate clusters for various window size
    logger("=== Clustering ===")
    clusters, _windowSize = bam2clusters(bam, fasta, outdir, windowSize, mapq, dpi, upto, verbose)

    # use preselected window size
    if len(windowSize)==1 and windowSize[0]*1000 in _windowSize:
        selected = _windowSize.index(windowSize[0]*1000)
    else:
        selected = np.argmin(map(len, clusters)) # check also number of contigs & cummulative size
        
    # choose window with least scaffolds
    logger("=== Scaffolding ===")
    _clusters, _windowSize = clusters[selected], _windowSize[selected]
    logger("Selected %s clusters from window %sk"%(len(_clusters), _windowSize/1000))
    #return
    
    for _windowSize2 in windowSize2:
        logger("=== %sk windows ==="%_windowSize2)
        logger("Constructing scaffolds...")
        scaffolds = get_scaffolds(bam, fasta.name, _clusters, mapq, minWindows, threads, verbose, _windowSize2)

        logger("Reporting %s scaffolds..."%len(scaffolds))
        outbase = outdir+"/%sk.%sk"%(_windowSize/1000, _windowSize2)
        fastafn = report_scaffolds(outbase, scaffolds, faidx)
    
    logger("Done!")

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--version', action='version', version='1.01b')   
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose")    
    parser.add_argument("-i", "--bam", nargs="+", help="BAM file(s)")
    parser.add_argument("-f", "--fasta", type=file, help="Genome FastA file")
    parser.add_argument("-o", "--outdir", required=1, help="output name")
    parser.add_argument("-w", "--windowSize", nargs="+", default=[1000, 500, 100, 50], type=int,
                        help="window size in kb used for karyotyping [%(default)s]")
    parser.add_argument("-z", "--windowSize2", default=[5, 10, 50], type=int, nargs="+", 
                        help="window size in kb used for scaffolding [%(default)s]")
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
    