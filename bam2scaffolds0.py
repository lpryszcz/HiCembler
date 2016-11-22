#!/usr/bin/env python
desc="""Report scaffolds by joining contigs based on contact matrix from BAM file. 

TBD:
- allow for fitting contig into gap within large contig?
- distance estimation?
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
from array2scaffolds import transform, normalize_rows, get_name, get_clusters, array2tree, tree2scaffold, _average_reduce_2d, report_scaffolds, getNewick
from FastaIndex import FastaIndex

# update sys.path & environmental PATH
root = os.path.dirname(os.path.abspath(sys.argv[0]))
src = ["bin", "bin/snap", "bin/sinkhorn_knopp"]
paths = [os.path.join(root, p) for p in src]
sys.path = paths + sys.path
os.environ["PATH"] = "%s:%s"%(':'.join(paths), os.environ["PATH"])

from sinkhorn_knopp import sinkhorn_knopp

def normalize(d):
    """Return fully balanced matrix"""
    sk = sinkhorn_knopp.SinkhornKnopp() #max_iter=100000, epsilon=0.00001)
    # make symmetric & normalise
    d += d.T
    d -= np.diag(d.diagonal()/2)
    d += 1
    d = sk.fit(d)
    return d

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

def get_longest(t, maxdist=6):
    """Return node having longest branch"""
    #n = sorted(t.traverse(), key=lambda n: 2*n.dist-t.get_distance(n), reverse=1)[0]
    n = t
    bestdist = 2*n.dist-n.get_distance(t)
    for _n in t.traverse():
        if _n.get_distance(t, topology_only=1) > maxdist:
            break
        if 2*_n.dist-_n.get_distance(t) > bestdist:
            n = _n
            bestdist = 2*_n.dist-_n.get_distance(t)
    return n
            
def get_subtrees(d, bin_chr, bin_position, method="ward", nchrom=1000, distfrac=0.75):
    maxtdist = 0
    i = 0
    subtrees = []
    for i in range(1, nchrom):
        names = ["%s %s"%(c, s) for c, (s, e) in zip(bin_chr, bin_position)]
        #Z = sch.linkage(d[np.triu_indices(d.shape[0], 1)], method=method)
        #tree = sch.to_tree(Z, False)
        #t = ete3.Tree(getNewick(tree, "", tree.dist, names))
        t = array2tree(d, names, method=method)
        
        tname, tdist = t.get_farthest_leaf()#[1]
        if maxtdist < tdist:
            maxtdist = t.get_farthest_leaf()[1]
        # break if small subtree
        if tdist < maxtdist*distfrac:
            break
          
        # get longest branch
        n = get_longest(t)
        pruned = n.get_leaf_names()         
        subtrees.append(pruned)
        # prune array        
        indices = [_i for _i, name in enumerate(names) if name not in set(pruned)]
        d = d[indices, :]
        d = d[:, indices]
        bin_chr = bin_chr[indices]
        bin_position = bin_position[indices, :]
            
    if i:    
        subtrees.append(t.get_leaf_names())
    return subtrees
    
def get_clusters(outbase, d, contig2size, bin_chr, bin_position, method="ward", nchrom=1000, distfrac=0.75):
    """Return clusters for given distance matrix"""
    maxtdist = 0
    d = transform(d)
    subtrees = get_subtrees(d, bin_chr, bin_position)

    logger(" Assigning contigs to %s clusters..."%len(subtrees))
    total = correct = 0
    contig2cluster = {get_name(c): Counter() for c in np.unique(bin_chr)}
    for i, subtree in enumerate(subtrees, 1):
        c = Counter(map(get_name, subtree))
        total += len(subtree)
        correct += c.most_common(1)[0][1]
        # poplate contig2clustre
        for k, v in c.iteritems():
            if not k: continue
            contig2cluster[get_name(k)][i] += v
    logger("  %s / %s [%.2f%s]"%(correct, total, 100.*correct/total, '%'))

    logger(" Weak assignments...")
    clusters = [[] for i in range(len(subtree))]
    withoutCluster, weakCluster = [], []
    for c, counter in contig2cluster.iteritems():
        if not counter:
            withoutCluster.append(c)
            continue
        # get major cluster
        clusteri, count = counter.most_common(1)[0]
        mfrac = 1. * count / sum(counter.itervalues())
        clusters[clusteri].append(c)
        if mfrac<.66:
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
    logger("  %3s bp in %s clusters."%(totsize, len(clusters)))
    return clusters        
    
def bam2clusters(bam, fasta, outdir, windowSize, mapq, dpi, upto, verbose):
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
    #npy = np.load(os.path.join(outdir,"100k.npz")); arrays = [npy[npy.files[0]], ]
    
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

        logger("Sinkhorn-Knopp normalisation...")
        d = normalize(d)

        # generate missing handles
        bin_chr, bin_position = [], []
        for c, s, e in _windows:
            bin_chr.append(c)
            bin_position.append((s, e))
        bin_chr = np.array(bin_chr)
        bin_position = np.array(bin_position)

        logger("Assigning contigs to clusters/scaffolds...")
        _clusters = get_clusters(outfn, d, contig2size, bin_chr, bin_position)
        clusters.append(_clusters)
    return clusters, windowSize

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
    # make symmetric & normalise
    d += normalize(d)
    
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

def get_scaffolds(bam, fasta, clusters, mapq, minWindows, threads, verbose, windowSize=10):
    """Return scaffolds"""
    scaffolds = []
    faidx = FastaIndex(fasta)
    contig2size = {c: stats[0] for c, stats in faidx.id2stats.iteritems()}
    if threads>1:
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
    parser.add_argument('--version', action='version', version='1.0b')   
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
    