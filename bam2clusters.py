#!/usr/bin/env python
desc="""Cluster contigs based on contact matrix from BAM file.

TBD:
- better clustering
- improve speed
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
from FastaIndex import FastaIndex
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
#from bam2scaffolds import normalize_rows, normalize_diagional

transform = lambda x: np.log(np.max(x+1))-np.log(x+1)

# update sys.path & environmental PATH
root = os.path.dirname(os.path.abspath(sys.argv[0]))
src = ["bin", "bin/snap", "bin/sinkhorn_knopp"]
paths = [os.path.join(root, p) for p in src]
sys.path = paths + sys.path
os.environ["PATH"] = "%s:%s"%(':'.join(paths), os.environ["PATH"])

from sinkhorn_knopp import sinkhorn_knopp

def normalize(d, max_iter=1000, epsilon=0.00001):
    """Return fully balanced matrix"""
    sk = sinkhorn_knopp.SinkhornKnopp(max_iter=max_iter, epsilon=epsilon)
    # make symmetric & normalise
    d += d.T
    d -= np.diag(d.diagonal()/2)
    d += 1
    vmax = d.max()
    d = sk.fit(d/vmax)
    d *= vmax / d.max()
    return d

def normalize_diagional(d, bin_chr, bin_position):
    """Return symmetric and diagonal normalised matrix"""
    logger("Diagonal normalisation...")    
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
    n2 = np.mean(d.diagonal()) / d.diagonal()
    d *= (d*n2).T*n2 #np.mean(d.diagonal()) / d.diagonal()
    return d, bin_chr, bin_position
    
def logger(message, log=sys.stdout):
    """Log messages"""
    memory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    log.write("[%s] %s    [memory: %6i Mb]\n"%(datetime.ctime(datetime.now()), message, memory))
    
def contigs2windows(fasta, minSize=2000, verbose=0):
    """Return one window per contig, filtering out contigs shorter than minSize"""
    genomeSize = 0
    windows, skipped, chr2window, contig2size = [], [], {}, {}
    faidx = FastaIndex(fasta)
    for i, c in enumerate(faidx, 1):
        if i%1e5 == 1:
            sys.stderr.write(' %s   \r'%i)
        # add contig size
        size = faidx.id2stats[c][0]
        contig2size[c] = size    
        # skip short contigs    
        if size < minSize:
            skipped.append(size)
            continue
        # store chr2window info and append window
        chr2window[c] = len(windows)
        windows.append((c, 0, size))
        # update genomeSize
        genomeSize += size
    if verbose:
        logger(' %s bases in %s contigs divided in %s windows. '%(faidx.genomeSize, len(faidx), len(windows)))
        if skipped:
            logger('  %s bases in %s contigs skipped.'%(sum(skipped), len(skipped)))
    return windows, chr2window, contig2size
    
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

def _bam2array(args):
    """Parse pysam alignments and return matching windows"""
    _bam, contigs, mapq, contig2size, chr2window, minSize = args 
    data = Counter() #for ii in range(len(windowSize))]
    regions = []
    for c in contigs: 
        regions += ["%s:%s-%s"%(c, 1, minSize), "%s:%s-%s"%(c, contig2size[c]-minSize, contig2size[c])]
    proc = _get_samtools_proc(_bam, mapq, regions)
    # process reads from given library
    for line in proc.stdout:
        # unload data
        ref1, start1, mapq, cigar, ref2, start2, insertSize, seq = line.split('\t')[2:10]
        start1, start2, seqlen = int(start1), int(start2), len(seq)
        # update ref if alignment in the same chrom
        if ref2 == "=":
            ref2 = ref1
        # skip if contig not present in array # or not in contig end regions
        if ref1 not in chr2window or ref2 not in chr2window or \
           start2 > minSize and start2 < contig2size[ref2]-minSize-len(seq):
            continue
        # get windows
        w1, w2 = chr2window[ref1], chr2window[ref2]
        # matrix is symmetrix, so make sure to add only to one part
        if w2 < w1:
            w1, w2 = w2, w1
        # update contact array
        data[(w1, w2)] += 1
                
    # stop subprocess
    proc.terminate()
    return data

def bam2array_multi(windows, contig2size, chr2window, bam, mapq, upto=0,
                    regions=[], verbose=1, minSize=2000, threads=4):
    """Return contact matrix based on BAM"""
    arrays = np.zeros((len(windows), len(windows)), dtype="float32")# for w in windows]
    i = 1
    # init empty array
    for _bam in bam:
        if verbose:
            logger(" %s"%_bam)
        if not regions:
            regions = chr2window.keys(); 
        args = [(_bam, regions[ii::threads], mapq, contig2size, chr2window, minSize) for ii in range(threads)]
        p = Pool(threads)
        for wdata in p.imap_unordered(_bam2array, args):#, chunksize=100):
            for (w1, w2), c in wdata.iteritems():
                arrays[w1][w2] += c
                i += c
            #i = arrays[0].sum()
            if upto and i>upto:
                break
            sys.stderr.write(" %s \r"%i)
    if verbose:
        logger(" %s alignments parsed"%i)
    return arrays
    
def getNewick(node, newick, parentdist, leaf_names):
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = getNewick(node.get_left(), newick, node.dist, leaf_names)
        newick = getNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick
    
def distance_matrix2tree(Z, names):
    """Return tree representation for distance matrix"""
    n = Z.shape[0]+1
    i2n = [0] * (2*n - 1)
    t = ete3.Tree()
    for i, (idx1, idx2, dist, sample_count) in enumerate(Z):
        idx1, idx2 = int(idx1), int(idx2)
        # create Tree object for tips / leaves
        if idx1 < n:
            i2n[idx1] = ete3.Tree(name=names[idx1])
        if idx2 < n:
            i2n[idx2] = ete3.Tree(name=names[idx2])
        # create new node
        t = ete3.Tree()
        # normalise distance
        dist1 = dist - i2n[idx1].get_farthest_leaf()[1]
        dist2 = dist - i2n[idx2].get_farthest_leaf()[1]
        # add children
        t.add_child(i2n[idx1], dist=dist1)
        t.add_child(i2n[idx2], dist=dist2)
        # store
        i2n[n + i] = t
    return t

def mylayout(node):
    # don't show circles
    node.img_style["size"] = 0
    node.img_style["vt_line_width"] = 0    
    # If node is a leaf, show aligned node name
    if node.is_leaf():
        nameFace = ete3.faces.TextFace(node.name, fsize=8)
        ete3.faces.add_face_to_node(nameFace, node, column=0, aligned=True)
        
def array2tree(d, names, outbase="", method="ward"):
    """Return tree representation for array"""
    '''# cluster
    Z = sch.linkage(d[np.triu_indices(d.shape[0], 1)], method=method)
    t = distance_matrix2tree(Z, names); print len(set(t.get_leaf_names())), len(names)
    '''
    Z = fastcluster.linkage(d[np.triu_indices(d.shape[0], 1)], method=method)
    tree = sch.to_tree(Z, False)
    t = ete3.Tree(getNewick(tree, "", tree.dist, names))
    #'''
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

def get_name(contig):
    return contig.split()[0].split('|')[-1].split('.')[0]
    
def report_tree(nw, fasta, mind=3, maxd=5):
    """Save tree pdf"""
    # load contig2chrom
    c2chr = {}
    if os.path.isfile("%s.bed"%fasta):
        c2chr = {l.split('\t')[3]: get_name(l.split('\t')[0]) for l in open("%s.bed"%fasta)}
    # load tree
    t = ete3.Tree(nw)
    # truncate branches
    for i, n in enumerate(t.traverse(), 1):
        dist = t.get_distance(n, topology_only=1)
        chrs = Counter(c2chr[c] if c in c2chr else c for c in n.get_leaf_names())
        if dist > mind and len(chrs)==1 or maxd and dist>maxd:
            n.leaves = n.get_leaf_names()
            n.chrs = chrs
            n.name="%s %s chrs %s leaves"%(chrs.most_common(1)[0][0], len(chrs), len(n))
            for _n in n.get_children():
                n.remove_child(_n)

    ts = ete3.TreeStyle()
    ts.show_leaf_name = False
    ts.layout_fn = mylayout
    t.render('truncated.tree.pdf', tree_style=ts)
    
def get_subtrees(d, bin_chr, bin_position, method="ward", nchrom=1000, distfrac=0.4, fasta=""):
    """Return contings clustered into scaffolds
    fastcluster is slightly faster than scipy.cluster.hierarchy and solve the issue: http://stackoverflow.com/a/40738977/632242
    """
    subtrees = [] 
    i = maxtdist = 0
    t = array2tree(d, bin_chr)#, "clustering.tree")
    report_tree(t.write(), fasta)
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

def bam2clusters(bam, fasta, outdir, minSize=2000, mapq=10, threads=4, dpi=100, upto=0, verbose=1, method="ward"):
    """Return clusters computed from from windowSizes"""
    logger("=== Clustering ==="); print minSize, threads
    outbase = os.path.join(outdir, "auto")
    
    # load clusters
    fname = outbase + ".clusters.tab_"
    if os.path.isfile(fname):
        clusters = [l[:-1].split('\t') for l in open(fname)]
        logger("  %s clusters loaded."%len(clusters))
        return clusters

    # get windows
    windows, chr2window, contig2size = contigs2windows(fasta, minSize, verbose)
        
    # generate missing handles
    bin_chr, bin_position = [], []
    for c, s, e in windows:
        bin_chr.append(c)
        bin_position.append((s, e))
    bin_chr = np.array(bin_chr)
    bin_position = np.array(bin_position)

    if not os.path.isfile(outbase+".npz"): #.balanced
        # get array from bam
        logger("Parsing BAM...")
        d = bam2array_multi(windows, contig2size, chr2window, bam,  mapq, upto, \
                            threads=threads, minSize=minSize)

        # save windows, array and plot
        logger("Saving array...")
        with gzip.open(outbase + ".windows.tab.gz", "w") as out:
            out.write("\n".join("\t".join(map(str, w)) for w in windows)+"\n")
        with open(outbase+".npz", "w") as out:
            np.savez_compressed(out, d)

        '''# make symmetric & normalise
        logger("Balancing array...")
        d = normalize(d)
        #d, bin_chr, bin_position = normalize_diagional(d, bin_chr, bin_position)
        #d, bin_chr, bin_position = normalize_window_size(d, bin_chr, bin_position, _windowSize)
        logger(" saving fully balanced array...")
        with open(outbase+".balanced.npz", "w") as out:
            np.savez_compressed(out, d)'''
    # load from file
    else:
        npy = np.load(outbase+".npz") #balanced.
        d = npy[npy.files[0]]
        
    logger("Balancing array...")
    d = normalize(d, 1)
    #d, bin_chr, bin_position = normalize_diagional(d, bin_chr, bin_position)
    #d, bin_chr, bin_position = normalize_window_size(d, bin_chr, bin_position, _windowSize)
    
    # get clusters on transformed matrix
    logger("Clustering..."); print d.max(), d.mean(), method
    clusters = get_subtrees(transform(d), bin_chr, bin_position, method, fasta=fasta)

    # skip empty clusters
    clusters = filter(lambda x: x, clusters)
    totsize = contigs = 0
    outfn = outbase+".clusters.tab"
    logger("Reporting %s clusters to %s ..."%(len(clusters), outfn))
    with open(outfn, "w") as out:
        for i, cluster in enumerate(clusters, 1):
            clSize = sum(contig2size[c] for c in cluster)
            totsize += clSize
            contigs += len(cluster)
            out.write("\t".join(cluster)+"\n")
    logger("  %3s bp in %s clusters generated from %s contigs."%(totsize, len(clusters), contigs))
    
    return clusters

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--version', action='version', version='1.01d')   
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose")    
    parser.add_argument("-i", "--bam", nargs="+", help="BAM file(s)")
    parser.add_argument("-f", "--fasta", type=file, help="Genome FastA file")
    parser.add_argument("-o", "--outdir", required=1, help="output name")
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
    bam2clusters(o.bam, o.fasta, o.outdir, o.minSize, o.mapq, o.threads, o.dpi, o.upto, o.verbose)
        
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
    