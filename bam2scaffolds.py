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

import ete3, glob, gzip, os, pysam, resource, subprocess, sys
import scipy.cluster.hierarchy as sch
import fastcluster
import numpy as np
from collections import Counter
from datetime import datetime
from multiprocessing import Pool
from fastq2array import logger, fasta2windows, get_window, plot
from array2scaffolds import transform, normalize_rows, get_name, get_clusters, array2tree, _average_reduce_2d, report_scaffolds, getNewick
from FastaIndex import FastaIndex
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

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

def normalize_window_size(d, bin_chr, bin_position):
    """Return symmetric and normalised matrix by window size"""
    # normalise by window size
    sizes = np.diff(bin_position, axis=1)[:, 0]
    c = Counter(sizes)
    windowSize, occurencies = c.most_common(1)[0]#; print windowSize
    d *= 1. * windowSize / sizes    
    return d, bin_chr, bin_position    

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

def bam2array(windows, windowSize, chr2window, bam, mapq, upto=1e7, regions=[], verbose=1, threads=1):
    """Return contact matrix based on BAM"""
    arrays = [np.zeros((len(w), len(w)), dtype="float32") for w in windows]
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
    
def _bam2array(args):
    """Parse pysam alignments and return matching windows"""
    _bam, regions, mapq, windowSize, chr2window = args 
    data = [Counter() for ii in range(len(windowSize))]
    proc = _get_samtools_proc(_bam, mapq, regions)
    # process reads from given library
    for line in proc.stdout:
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
                data[ii][(w1, w2)] += 1
                
    # stop subprocess
    proc.terminate()
    return data

def bam2array_multi(windows, windowSize, chr2window, bam, mapq, upto=1e7, regions=[], verbose=1, threads=1):
    """Return contact matrix based on BAM"""
    arrays = [np.zeros((len(w), len(w)), dtype="float32") for w in windows]    
    i = 1
    # init empty array
    for _bam in bam:
        if verbose:
            logger(" %s"%_bam)
        if not regions:
            regions = chr2window[0].keys()#; print regions[:10]
        args = [(_bam, [region], mapq, windowSize, chr2window) for region in regions]
        p = Pool(threads)
        for wdata in p.imap_unordered(_bam2array, args, chunksize=100):
            for ii, wdataii in enumerate(wdata):
                #arrays[ii] += wdataii
                for (w1, w2), c in wdataii.iteritems():
                    arrays[ii][w1][w2] += c
                    if not ii: i += c
            #i = arrays[0].sum()
            if upto and i>upto:
                break
            sys.stderr.write(" %s \r"%i)
    if verbose:
        logger(" %s alignments parsed"%i)
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
            clSize = sum(contig2size[c] for c in cluster)
            totsize += clSize
            out.write("\t".join(cluster)+"\n")
    logger("  %3s bp in %s clusters generated from %s contigs."%(totsize, len(clusters), len(contig2cluster)))
    return clusters        
    
def bam2clusters(bam, fasta, outdir, windowSize, mapq, threads, dpi, upto, verbose):
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
    arrays = bam2array_multi(windows, windowSize, chr2window, bam,  mapq, upto, threads=threads)
    
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
        #d, bin_chr, bin_position = normalize_diagional(d, bin_chr, bin_position)
        d, bin_chr, bin_position = normalize_window_size(d, bin_chr, bin_position)

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
    
def tree2scaffold(t, d, bin_chr, bin_position, contig2indices, minWindows):
    """Scaffold contigs based on distance matrix and tree."""
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

def contigs2scaffold(args):
    """Combine contigs into scaffold"""
    i, bam, fasta, windowSize, contigs, mapq, minWindows, params = args
    
    # get windows
    data = fasta2windows(fasta, windowSize, verbose=0, skipShorter=1, contigs=set(contigs), filterwindows=0)
    windowSize, windows, chr2window, base2chr, genomeSize = data

    # get array from bam
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
    
    logger(" cluster_%s with %s windows for %s out of %s contigs"%(i, d.shape[0], len(np.unique(bin_chr)), len(contigs)))

    # skip if empty array or report without scaffolding if only one contig
    if not d.shape[0]:
        return []
    if len(np.unique(bin_chr)) == 1:
        return [(bin_chr[0], 0)]
        
    # make symmetric & normalise
    d += d.T
    d -= np.diag(d.diagonal()/2)
    d = normalize_rows(d)

    # get tree on reduced matrix
    t = array2tree(transform(_average_reduce_2d(d, bin_chr)), np.unique(bin_chr))
    
    # get contig with indices
    contig2indices = get_contig2indices(bin_chr)
    
    # get scaffold
    scaffold = tree2scaffold(t, d, bin_chr, bin_position, contig2indices, minWindows)

    # estimate & add gaps
    faidx = FastaIndex(fasta)
    contig2size = {c: stats[0] for c, stats in faidx.id2stats.iteritems()}
    scaffold = scaffold2gaps(d, bin_chr, bin_position, scaffold, params, contig2indices, contig2size)
    return scaffold

def scaffold2gaps(d, bin_chr, bin_position, scaffold, params, contig2indices, contig2size):
    """Add gaps to scaffolds"""
    d = normalize(d)
    scaffold_with_gaps = []
    for i, (c, o) in enumerate(scaffold[:-1]):
        c2, o2 = scaffold[i+1]
        contacts = d[contig2indices[c]][:, contig2indices[c2]].reshape(len(contig2indices[c])*len(contig2indices[c2]))
        distances = [distance_func(_c, *params) for _c in contacts] # skip outliers ie 0.9 percentile
        gapSize = int(np.median(distances)*1000 - contig2size[c]/2 - contig2size[c2]/2) 
        print gapSize, contig2size[c], contig2size[c2], np.mean(distances), np.std(distances), distances
        scaffold_with_gaps.append((c, o, gapSize))
    # add last contig without gap
    c, o = scaffold[-1]
    scaffold_with_gaps.append((c, o, 0))
    return scaffold_with_gaps

def distance_func(x, a, b):
    """Function to calculating distance from normalised contact frequency"""
    return (1./(a*x)) **(1./b)
    
def contact_func(x, a, b):
    """Function to calculating normalised contact frequency from distance"""
    return 1./(a * x ** b)

def estimate_distance_parameters(out, bam, mapq, contig2size, windowSize=10, skipfirst=5, icontigs=25, upto=1e5):
    """Return estimated parameters"""
    logger(" Estimating distance parameters...")
    i = 0
    d2c = {0: []}
    longest_contigs = sorted(contig2size, key=lambda x: contig2size[x], reverse=1)
    maxdist = 0.75 * contig2size[longest_contigs[0]] / 1000.
    for c in longest_contigs[:icontigs]:
        s = contig2size[c]
        sys.stderr.write(' %s %s bp   \r'%(c, s))
        n = s / windowSize + 1 # since windowSize is in kb
        positions = range(windowSize, n*windowSize, windowSize)
        arrays = [np.zeros((n, n), dtype='float32')]
        chr2window = [{c: 0}]
        arrays = bam2array(arrays, [windowSize], chr2window, bam, mapq, regions=[c], upto=upto, verbose=0)
        d = arrays[0]
        i += d.sum()
        d += d.T
        d -= np.diag(d.diagonal()/2)
        d = normalize(d) #normalize_rows(d)
        
        d2c[0] += [d[i][i] for i in range(d.shape[0]-1)]
        for i in range(len(positions)-1):
            for j in range(i, len(positions)):
                dist = (positions[j]-positions[i]) / 1000 
                if dist not in d2c:
                    d2c[dist] = []
                d2c[dist].append(d[i][j])
        if i>upto:
            break

    # plot up to 0.75x max dist
    dists = np.array([d for d in sorted(d2c)[:30] if d<maxdist])
    contacts = [d2c[d] for d in dists]

    plt.title("Normalised HiC contacts vs distance")
    plt.boxplot(contacts, 1, '', positions=dists, widths=.75*dists[1])
    plt.xticks(rotation=90)

    plt.xlabel("Genomic distance [kb]")
    plt.ylabel("Normalised contacts")
    plt.xlim(xmin=-dists[1])
    plt.ylim(ymin=0)

    x = dists
    yn = np.array([np.median(c) for c in contacts])
    step = dists[1]/4.
    xs = np.arange(0, max(x)+step, step)

    params, pcov = curve_fit(contact_func, x[skipfirst:], yn[skipfirst:])
    plt.plot(xs[1:], contact_func(xs[1:], *params), 'r-', label="Fit\n$ y = 1 / {%0.2f x^ {%0.2f}} $"%tuple(params))

    outfn1, outfn2 = out+".distance_fit.png", out+".distance_fit.zoom.png"
    plt.legend(fancybox=True, shadow=True)
    plt.savefig(outfn1)
    plt.ylim(ymin=0, ymax=yn[0]/30.)
    plt.savefig(outfn2)
    logger("  fit stored as %s"%outfn1)
    return params
    
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
    
def bam2scaffolds(bam, fasta, outdir, windowSize, windowSize2, mapq, threads, dpi, upto, minWindows, verbose):
    """Report scaffolds based on BAM file"""
    faidx = FastaIndex(fasta)

    # calculate clusters for various window size
    logger("=== Clustering ===")
    clusters, _windowSize = bam2clusters(bam, fasta, outdir, windowSize, mapq, threads, dpi, upto, verbose)

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
        outbase = outdir+"/%sk.%sk"%(_windowSize/1000, _windowSize2)
        
        logger("=== %sk windows ==="%_windowSize2)
        logger("Constructing scaffolds...")
        scaffolds = get_scaffolds(outbase, bam, fasta.name, _clusters, mapq, minWindows, threads, verbose, _windowSize2)

        logger("Reporting %s scaffolds..."%len(scaffolds))
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
    parser.add_argument("-w", "--windowSize", nargs="+", default=[100, 50, 20, 10], type=int,
                        help="window size in kb used for karyotyping [%(default)s]")
    parser.add_argument("-z", "--windowSize2", default=[5, 10, 50, 2], type=int, nargs="+", 
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
    