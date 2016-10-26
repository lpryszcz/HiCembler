#!/usr/bin/env python
desc="""Report scaffolds by joining contigs based on contact matrix. 

TBD:
- estimate distances between contigs in scaffold.
- multithread scaffolding
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Bratislava, 25/10/2016
"""

import ete3, gzip, os, resource, sys
import scipy.cluster.hierarchy as sch
import numpy as np
from collections import Counter
from datetime import datetime
from fastq2array import logger
from FastaIndex import FastaIndex

def normalize_rows(a):
    """Normalise rows of give array."""
    rows, cols = a.shape
    maxv = a.sum(axis=0).max()
    for i in xrange(rows):
        # only if any signal
        if a[i].max():
            a[i] *= maxv/a[i].sum()
    return a
    
def load_matrix(fname, chrs=[], remove_shorter=False, verbose=0):
    """Load Hi-C interaction matrix from numpy dump
    generated by fastq2array.py. 
     
    Returns:
    d: data matrix over the selected set of chromosomes.
    bin_chr: list of chr index assignment of each bin.
    bin_position: start and end position of each bin
    """
    # load array
    npy = np.load(fname)
    d = npy[npy.files[0]]
    
    # load windows
    windowfn = fname[:-4]+'.windows.tab.gz' 
    bin_chr = []
    bin_position = [] 
    for i, l in enumerate(gzip.open(windowfn)):
        chrom, start, end = l[:-1].split('\t')
        bin_chr.append(chrom)
        bin_position.append(map(int, (start, end)))
        
    # chromosome array
    bin_position = np.array(bin_position)
    bin_chr = np.array(bin_chr)

    # calculate genome size and contig2size
    contig2size = {get_name(c): 0 for c in np.unique(bin_chr)}
    for c, (s, e) in zip(bin_chr, bin_position):
        contig2size[get_name(c)] += e-s

    #''' # eliminate 
    if remove_shorter:
        c = Counter(np.diff(bin_position, axis=1)[:, 0])
        windowSize, occurencies = c.most_common(1)[0]
        if verbose:
            sys.stderr.write(" most common window: %s bp [%5.2f%s]\n"%(windowSize, occurencies*100./len(bin_chr), '%'))
        valid_rowcols = ~(np.diff(bin_position, axis=1)[:, 0]!=windowSize)
        d = d[valid_rowcols, :]
        d = d[:, valid_rowcols]
        bin_chr = bin_chr[valid_rowcols]
        bin_position = bin_position[valid_rowcols, :]
    '''
    # contigs shorter than 3x window size    
    if remove_shorter and not chrs:
        c = Counter(np.diff(bin_position, axis=1)[:, 0])
        windowSize, occurencies = c.most_common(1)[0]
        # get contigs having size 3x bigger than windowSize
        contig2size = {c: s for c, s in filter(lambda x: x[1]>2*windowSize, contig2size.iteritems())}
        chrs = contig2size.keys()
        if verbose:
            sys.stderr.write(" keeping %s contigs having size > 2*%s ...\n"%(len(chrs), windowSize))
    #'''    
    # keep only relevant chromosomes
    if chrs:
        indices = np.any(bin_chr[None].T == chrs, 1)
        d = d[indices, :]
        d = d[:, indices]
        bin_chr = bin_chr[indices]
        bin_position = bin_position[indices, :]
        
    # make symmetric at the end to save memory and time
    d += d.T
    d -= np.diag(d.diagonal()/2)
    d = normalize_rows(d)
        
    return d, bin_chr, bin_position, contig2size

def distance_matrix2tree(Z, names):
    """Return tree representation for distance matrix"""
    n = Z.shape[0]+1
    i2n = {}
    idx = 0
    t = ete3.Tree()
    for i, (idx1, idx2, dist, sample_count) in enumerate(Z, 1):
        idx1, idx2 = int(idx1), int(idx2)
        # create Tree object for tips / leaves
        if idx1 < n:
            i2n[idx1] = ete3.Tree(name=names[idx1], dist=0)
        if idx2 < n:
            i2n[idx2] = ete3.Tree(name=names[idx2], dist=0)
        # create new node
        t = ete3.Tree(dist=0)
        # normalise distance
        dist -= max(i2n[idx1].get_farthest_leaf()[1], i2n[idx2].get_farthest_leaf()[1])
        # add children
        t.add_child(i2n[idx1], dist=dist)
        t.add_child(i2n[idx2], dist=dist)
        # store
        i2n[n+idx] = t
        idx += 1
    return t
        
def array2tree(d, names, outbase="", method="ward"):
    """Return tree representation for array"""
    # cluster
    Z = sch.linkage(d[np.triu_indices(d.shape[0], 1)], method=method)
    
    # get ete Tree
    t = distance_matrix2tree(Z, names)
    
    # save tree & newick
    if outbase:
        pdf, nw = outbase+".nw.pdf", outbase+".nw"
        with open(nw, "w") as out:
            out.write(t.write())
            
        ts = ete3.TreeStyle()
        ts.show_leaf_name = False
        ts.layout_fn = mylayout
        t.render(pdf, tree_style=ts)
        
    return t

def get_name(contig):
    return contig.split()[0] 

def mylayout(node):
    # don't show circles
    node.img_style["size"] = 0
    node.img_style["vt_line_width"] = 0    
    # If node is a leaf, show aligned node name
    if node.is_leaf():
        nameFace = ete3.faces.TextFace(node.name, fsize=8)
        ete3.faces.add_face_to_node(nameFace, node, column=0, aligned=True)

def get_shuffled(d, bin_chr, bin_position, seed=0):
    prng = np.random.RandomState(seed=seed)
    perm = prng.permutation(d.shape[0])
    inv_perm = np.argsort(perm)
    return d[perm, :][:, perm], bin_chr[perm], bin_position[perm]            

def get_clusters(outbase, infile, t, contig2size, bin_chr):
    """Return clusters from tree"""
    # generate clusters
    subtrees=[]
    while len(t)>2:
        n, dist = t.get_farthest_leaf()
        dists = [a.dist for a in n.get_ancestors()]
        # get ancestor with the longest branch length
        ai = dists.index(max(dists))
        a = n.get_ancestors()[ai]
        if n.name:
            c = Counter(_n.name.split('.')[0] for _n in a)
            subtrees.append(a)
        p = n.get_ancestors()[ai+1]
        p.remove_child(a)

    logger("Assigning contigs to %s clusters..."%len(subtrees))
    total = correct = 0
    contig2cluster = {get_name(c): Counter() for c in np.unique(bin_chr)}
    for i, subtree in enumerate(subtrees, 1):
        c = Counter(get_name(_n.name) for _n in subtree if _n.name)
        total += len(subtree)
        correct += c.most_common(1)[0][1]
        # poplate contig2clustre
        for k, v in c.iteritems():
            if not k: continue
            contig2cluster[get_name(k)][i] += v
    print " %s / %s [%.2f%s]"%(correct, total, 100.*correct/total, '%')

    #logger("Weak assignments...")
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
    print "  %s bp in %s contigs without assignment."%(sum(contig2size[c] for c in withoutCluster), len(withoutCluster))
    print "  %s bp in %s contigs having weak assignment."%(sum(contig2size[c] for c in weakCluster), len(weakCluster))

    clusters = filter(lambda x: x, clusters)
    outfile = outbase+".clusters.tab"
    logger("Reporting clusters to: %s ..."%outfile)
    totsize = 0
    with open(outfile, "w") as out:
        for i, cluster in enumerate(clusters, 1):
            clSize = sum(contig2size[c] for c in cluster)
            totsize += clSize
            print " %s %s bp in %s contigs" % (i, clSize, len(cluster))
            out.write("\t".join(cluster)+"\n")
    print "%3s bp in %s clusters!"%(totsize, len(clusters))
    return clusters

def get_reversed(scaffold):
    """Return reversed scaffold, updating orientation"""
    return [(name, not orientation) for name, orientation in reversed(scaffold)]

def get_indices(scaffold, contig2indices):
    """Return list with indices representing given scaffold"""
    indices = []
    for name, reverse in scaffold:
        if reverse:
            indices += reversed(contig2indices[name])
        else:
            indices += contig2indices[name]
    return indices
    
def join_scaffolds(scaffold1, scaffold2, d, contig2indices, minWindows=3):
    """Join two adjacent scaffolds"""
    indices1, indices2 = get_indices(scaffold1, contig2indices), get_indices(scaffold2, contig2indices)
    # skip contigs with less windows than minWindows
    (scaffold1, indices1), (scaffold2, indices2) = sorted(((scaffold1, indices1), (scaffold2, indices2)), key=lambda x: len(x[1]), reverse=1)
    if len(indices2) < minWindows:
        return scaffold1
    # get subset of array for scaffold1 and scaffold2
    _d = d[:, indices1+indices2][indices1+indices2, :]
    # get orientation: 0: s-s; 1: s-e; 2: e-s; 3: e-e
    ## contact values for start and end of two contigs are compared
    ## and the max value is taken as true contact
    ## this may be replaced by some ML function, as comparing only 1 window may be misleading!
    n1 = len(indices1)
    orientation = np.argmax(_d[:, [n1, -1]][[0, n1-1],:])
    # s - s
    if   orientation == 0:
        scaffold = get_reversed(scaffold1) + scaffold2
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
        
def tree2scaffold(t, d, bin_chr, bin_position, minWindows):
    """Scaffold contigs based on distance matrix and tree."""
    # get contig with indices
    contig2indices = {c: [] for c in np.unique(bin_chr)}
    for i, c in enumerate(bin_chr):
        contig2indices[c].append(i)
        
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
        #print n.scaffold

    # estimate distances
    
        
    return t.scaffold

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

def _contigs2scaffold(args):
    """Combine contigs into scaffold"""
    infile, contigs = args
    # get part of matrix for particular scaffold
    _d, _bin_chr, _bin_position, contig2size = load_matrix(infile, chrs=contigs, remove_shorter=True)
    sys.stderr.write(" %s with %s windows in %s contigs\r"%(i, _d.shape[0], len(contigs)))
    # get tree on reduced matrix
    t = array2tree(transform(_average_reduce_2d(_d, _bin_chr)), np.unique(_bin_chr)) 
    # get scaffold
    scaffold = tree2scaffold(t, _d, _bin_chr, _bin_position, minWindows)
    return scaffold

def clusters2scaffolds(clusters, infile, minWindows, transform):
    """Process clusters into scaffolds."""
    # easy for multithreading!
    scaffolds = []
    for i, contigs in enumerate(clusters, 1):
        # get part of matrix for particular scaffold
        _d, _bin_chr, _bin_position, contig2size = load_matrix(infile, chrs=contigs, remove_shorter=True)
        sys.stderr.write(" %s with %s windows in %s contigs\r"%(i, _d.shape[0], len(contigs)))
        # get tree on reduced matrix
        t = array2tree(transform(_average_reduce_2d(_d, _bin_chr)), np.unique(_bin_chr)) 
        # get scaffold
        scaffold = tree2scaffold(t, _d, _bin_chr, _bin_position, minWindows)
        scaffolds.append(scaffold)
    return scaffolds
    
def clusters2scaffolds1(clusters, infile, minWindows, transform):
    """Process clusters into scaffolds."""
    # easy for multithreading!
    scaffolds = []
    for i, contigs in enumerate(clusters, 1):
        # get part of matrix for particular scaffold
        _d, _bin_chr, _bin_position, contig2size = load_matrix(infile, chrs=contigs, remove_shorter=True)
        sys.stderr.write(" %s with %s windows in %s contigs\r"%(i, _d.shape[0], len(contigs)))
        # get tree on reduced matrix
        t = array2tree(transform(_average_reduce_2d(_d, _bin_chr)), np.unique(_bin_chr)) 
        # get scaffold
        scaffold = tree2scaffold(t, _d, _bin_chr, _bin_position, minWindows)
        scaffolds.append(scaffold)
    return scaffolds
            
def report_scaffolds(outbase, infile, scaffolds, faidx, w=60):
    """Save scaffolds"""
    totsize = 0
    fastafn = outbase+".scaffolds.fa"
    scaffoldfn = outbase+".scaffolds.tab"
    with open(fastafn, "w") as out, open(scaffoldfn, "w") as outscaffolds:
        for i, scaffold in enumerate(scaffolds, 1):
            seqs = []
            elements = []
            for c, o in scaffold:
                seqs.append(faidx.get_sequence(c, reverse=o))
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
    
def array2scaffolds(outbase, infile, infile2, fasta, minWindows, scaffolds=[]):
    """Return scaffolds computed for given matrix"""
    logger("Loading FastA...")
    faidx = FastaIndex(fasta)
    logger(" %s bp in %s contigs"%(faidx.genomeSize, len(faidx)))
    
    logger("Loading matrix from %s ..."%infile)
    d, bin_chr, bin_position, contig2size = load_matrix(infile, scaffolds=scaffolds, remove_shorter=True, verbose=1)
    logger(" matrix of %s windows for %s contigs summing %s bp"%(d.shape[0], len(contig2size), sum(contig2size.values())))

    # make sure all contigs from matrix are present in FastA
    diff = set(contig2size.keys()).difference(faidx)
    if diff:
        sys.stderr.write("[ERROR] %s / %s contigs are missing from provided FastA!\n"%(len(diff), len(contig2size)))
        sys.exit(1)

    #if shuffle: # don't use as positions assumed not shuffled matrix
    #    d, bin_chr, bin_position = get_shuffled(d, bin_chr, bin_position)

    logger("Calculating linkage matrix & tree...")
    transform = lambda x: np.log(np.max(x+1))-np.log(x+1)
    names = ["%s %7sk"%(get_name(c), s/1000) for c, (s, e) in zip(bin_chr, bin_position)]
    d = transform(d)
    t = array2tree(d, names, outbase)
    del d

    logger("Assigning contigs to clusters/scaffolds...") 
    clusters = get_clusters(outbase, infile, t, contig2size, bin_chr)

    logger("Constructing scaffolds...")
    scaffolds = clusters2scaffolds(clusters, infile2, minWindows, transform)

    logger("Reporting %s scaffolds..."%len(clusters))
    report_scaffolds(outbase, infile, scaffolds, faidx)
    
    return scaffolds

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--version', action='version', version='1.0b')   
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose")    
    parser.add_argument("-i", "--infile", required=True,
                        help="Contact matrix (.npz) used for chromosome assignment")
    parser.add_argument("-j", "--infile2", default='', 
                        help="Contact matrix (.npz) used for scaffolding [-i/--infile]")
    parser.add_argument("-o", "--outbase", default='', 
                        help="Output file base name [-i/--infile]")
    parser.add_argument("-f", "--fasta", required=True, type=file,
                        help="Contigs FastA file")
    parser.add_argument("-m", "--minWindows", default=3, type=int,
                        help="minimum number of windows per contigs [%(default)s]")

    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    infile = infile2 = outfile = o.infile
    if o.infile2:
        infile2 = o.infile2
    if o.outbase:
        outbase = o.outbase

    # first iteration
    logger("===== 1st ITERATION =====")
    scaffolds = array2scaffolds(outbase, infile, infile2, o.fasta, o.minWindows)

    # second iteration
    logger("===== 2nd ITERATION =====")
    scaffolds = array2scaffolds(outbase, infile, infile2, o.fasta, o.minWindows, scaffolds)    
    
if __name__=="__main__":
    t0 = datetime.now()
    main()
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
