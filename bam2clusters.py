#!/usr/bin/env python
desc="""Cluster contigs based on contact matrix from BAM file.

TBD:
- better clustering
- improve speed
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Bratislava, 27/10/2016
"""

# Force matplotlib to not use any Xwindows backend.
import matplotlib; matplotlib.use('Agg')

import glob, gzip, os, pickle, resource, subprocess, sys
import scipy.cluster.hierarchy as sch
import fastcluster
import numpy as np
from collections import Counter
from datetime import datetime
from multiprocessing import Pool
from FastaIndex import FastaIndex
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#transform = lambda x: np.log(np.max(x+1))-np.log(x+1)

# update sys.path & environmental PATH
root = os.path.dirname(os.path.abspath(sys.argv[0]))
src = ["bin", "bin/FastaIndex", "bin/idba/bin", "bin/snap", "bin/sinkhorn_knopp"]
paths = [os.path.join(root, p) for p in src]
sys.path = paths + sys.path
os.environ["PATH"] = "%s:%s"%(':'.join(paths), os.environ["PATH"])


def normalize(d, max_iter=1000, epsilon=0.00001):
    """Return fully balanced matrix"""
    from sinkhorn_knopp import sinkhorn_knopp
    sk = sinkhorn_knopp.SinkhornKnopp(max_iter=max_iter, epsilon=epsilon)
    # make symmetric & normalise
    d += d.T - np.diag(d.diagonal())
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
    
def contigs2windows(fasta, minSize=2000, verbose=0, genomeFrac=0.90):
    """Return one window per contig, filtering out contigs shorter than minSize"""
    genomeSize = 0
    windows, skipped, chr2window,  = [], [], {}
    # get N90 & update minSize if smaller
    faidx = FastaIndex(fasta)
    contig2size = {c: faidx.id2stats[c][0] for c in faidx}
    n90 = faidx.get_N_and_L(genomeFrac)
    if n90 > minSize:
        minSize = n90
        logger(" updated --minSize to N%s: %s"%(int(100*genomeFrac), minSize))
    # parse fasta & generate windows
    for i, c in enumerate(faidx, 1):
        if i%1e5 == 1:
            sys.stderr.write(' %s   \r'%i)
        # skip short contigs   
        size = contig2size[c]
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
    return minSize, windows, chr2window, contig2size
    
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

def get_window(chrom, start, length, windowSize, chr2window):
    """Return window alignment belongs to. """
    end = start + length
    # get middle position
    pos = int(start+round((end-start)/2))
    # get window
    window = pos / windowSize + chr2window[chrom]
    return window

def bam2arrays(windows, windowSize, chr2window, bam, mapq, upto=0, regions=[], verbose=1, threads=1):
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
    _bam, contigs, mapq, contig2size, chr2window, minSize = args 
    data = Counter() 
    regions = []
    for c in contigs: 
        regions += ["%s:%s-%s"%(c, 1, minSize), "%s:%s-%s"%(c, contig2size[c]-minSize, contig2size[c])]
    # process reads from given library
    c2dists = {c: [] for c in contigs}
    i = 0    
    proc = _get_samtools_proc(_bam, mapq, regions)
    for i, line in enumerate(proc.stdout, 1):
        # unload data
        ref1, start1, mapq, cigar, ref2, start2, isize, seq = line.split('\t')[2:10]
        start1, start2, seqlen = int(start1), int(start2), len(seq)
        # update ref if alignment in the same chrom
        if ref2 == "=":
            c2dists[ref1].append(int(isize))
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
    return data, i, c2dists

def bam2array(windows, contig2size, chr2window, bam, mapq=10, upto=0, \
              regions=[], minSize=2000, threads=4, verbose=1):
    """Return contact matrix and insert sizes from BAM
    Process only ends of contigs.
    """
    # init empty array
    array = np.zeros((len(windows), len(windows)), dtype="float32")
    if not regions: 
        regions = chr2window.keys()
    if threads > len(regions):
        threads = len(regions)
    # process regions / contigs
    i = 0
    c2dists = {}    
    p = Pool(threads)
    args = [(_bam, regions[ii::threads], mapq, contig2size, chr2window, minSize)
            for _bam in bam for ii in range(threads)]
    for wdata, algs, _c2dists in p.imap_unordered(_bam2array, args):
        i += algs
        for (w1, w2), c in wdata.iteritems():
            array[w1][w2] += c
        # this may be risky for very fragmented assemblies
        for c in _c2dists:
            if c not in c2dists:
                c2dists[c] = []
            c2dists[c] += _c2dists[c]
        if upto and i>upto:
            break
        sys.stderr.write(" %s \r"%i)
    if verbose:
        logger(" %s alignments parsed"%i)
    return array, c2dists

def distance_func(x, a, b, maxv=float('inf')):
    """Function to calculating distance from normalised contact frequency"""
    #if not x: return maxv
    return (1./(a*x)) **(1./b)

def contact_func(x, a, b):
    """Function to calculating normalised contact frequency from distance"""
    return 1./(a * x ** b)
    
def get_distances_contacts(bam, mapq, contig2size, windowSize, icontigs=5, upto=1e7):
    """Return estimated parameters."""
    logger(" Estimating distance parameters...")
    i = 0
    d2c = {0: []}
    longest_contigs = sorted(contig2size, key=lambda x: contig2size[x], reverse=1)
    maxdist = 0.75 * contig2size[longest_contigs[0]] / 1000.
    for c in longest_contigs[:icontigs]:
        s = contig2size[c]
        sys.stderr.write(' %s %s bp\t%s    \r'%(c, s, i))
        n = s / windowSize + 1 # since windowSize is in kb
        positions = range(windowSize, n*windowSize, windowSize)
        arrays = [np.zeros((n, n), dtype='float32')]
        chr2window = [{c: 0}]
        arrays = bam2arrays(arrays, [windowSize], chr2window, bam, mapq, regions=[c], upto=upto, verbose=0)
        d = arrays[0]
        i += d.sum()
        d += d.T - np.diag(d.diagonal())
        #d = normalize(d) #normalize_rows(d)
        
        d2c[0] += [d[i][i] for i in range(d.shape[0]-1)]
        for i in range(len(positions)-1):
            for j in range(i, len(positions)):
                dist = (positions[j]-positions[i]) / 1000 
                if dist not in d2c:
                    d2c[dist] = []
                d2c[dist].append(d[i][j])

    # plot up to 0.75x max dist
    dists = np.array([d for d in sorted(d2c) if d<maxdist]) #[:30]
    contacts = [d2c[d] for d in dists]
    return dists, contacts
    
def estimate_distance_parameters(outbase, bam=None, mapq=10, contig2size=None, windowSize=2000, \
                                 skipfirst=5, icontigs=5, c2dists={}, limit=30, upto=1e7):
    """Return estimated fit parameters.
    bam, contig2size & windowSize or dists & contacts have to be provided.
    """
    logger(" Estimating distance parameters...")
    # process bam
    if bam and not iSizes:
        dists, contacts = get_distances_contacts(bam, mapq, contig2size, windowSize, icontigs, upto)
    # or use pre-computed insert sizes from longest contigs
    else:
        dists = range(0, (limit+10)*windowSize, windowSize)
        contacts = [[] for i in range(len(dists)+10)]
        longest_contigs = sorted(contig2size, key=lambda x: contig2size[x], reverse=1)[:50]
        for c in longest_contigs: 
            for i, nc in  Counter(np.digitize(c2dists[c], dists)).iteritems():
                contacts[i-1].append(nc)
        # remove empty
        _dists, _contacts = dists[:-1], contacts[:-1]
        dists, contacts = [], []
        for d, c in zip(_dists, _contacts):
            if c:
                dists.append(d/1000.)
                contacts.append(c)
        
    # keep only smallest dists for estimation & plotting
    dists, contacts = dists[:limit], contacts[:limit]
    
    plt.figure()
    plt.title("HiC contacts vs distance")
    plt.boxplot(contacts, 1, '', positions=dists, widths=.75*dists[1])
    plt.xticks(rotation=90)

    plt.xlabel("Genomic distance [kb]")
    plt.ylabel("No. of contacts")
    plt.xlim(xmin=-dists[1])
    plt.ylim(ymin=0)
    plt.savefig(outbase+".distance.png")
    
    x = dists
    yn = np.array([np.median(c) for c in contacts])
    step = dists[1]/4.
    xs = np.arange(0, max(x)+step, step)

    params, pcov = curve_fit(contact_func, x[skipfirst:], yn[skipfirst:])#, maxfev=1000)
    outfn = outbase+".distance.params"
    with open(outfn, "w") as out:
        pickle.dump(params, out)
    # plot fit
    plt.plot(xs[1:], contact_func(xs[1:], *params), 'r-', label="Fit\n$ y = 1 / {%0.3f x^ {%0.3f}} $"%tuple(params))
    
    outfn1, outfn2 = outbase+".distance_fit.png", outbase+".distance_fit.zoom.png"
    plt.legend(fancybox=True, shadow=True)
    plt.savefig(outfn1)
    plt.ylim(ymin=0, ymax=yn[0]/30.)
    plt.savefig(outfn2)
    logger("  fit stored as %s"%outfn1)
    return params
    
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
    import ete3
    Z = fastcluster.linkage(d[np.triu_indices(d.shape[0], 1)], method=method)
    tree = sch.to_tree(Z, False)
    t = ete3.Tree(getNewick(tree, "", tree.dist, names))
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

''' 
def worker(args):
    d, prng, n, frac = args
    subset = int(round(n * frac))
    selection = prng.permutation(n)[:subset]
    iZ = fastcluster.linkage(d[selection, :][:, selection][np.triu_indices(subset, 1)], method="ward")
    return np.diff(iZ[:, 2])
'''

def cluster_contigs(outbase, d, bin_chr, bin_position, threads=4, frac=0.8,
                    iterations=50, maxnumchr=500, seed=0, dpi=100):
    """Estimate number of chromosomes and assign contigs to chromosomes"""
    logger(" Estimating number of chromosomes...")
    prng = np.random#.RandomState()
    Z = fastcluster.linkage(d[np.triu_indices(d.shape[0], 1)], method="ward")
    
    dZ = []
    n = d.shape[0]
    ''' 24 vs 31s on 4 cores - it's not worth as need to pickle large matrix and pass it through
    args = ((d, prng, n, frac) for i in range(iterations))
    p = Pool(threads)
    for i, _dZ in enumerate(p.imap_unordered(worker, args), 1):
        sys.stderr.write(' %s / %s  \r'%(i, iterations))
        dZ.append(_dZ)
    '''
    for i in range(iterations):
        sys.stderr.write(' %s / %s  \r'%(i+1, iterations))
        subset = int(round(n * frac))
        selection = prng.permutation(n)[:subset]
        iZ = fastcluster.linkage(d[selection, :][:, selection][np.triu_indices(subset, 1)], method="ward")
        dZ.append(np.diff(iZ[:, 2]))
    #'''
    # get number of chromosomes
    mean_step_len = np.mean(np.array(dZ), 0)
    nchr = np.argmax(mean_step_len[::-1][:maxnumchr]) + 2
    logger("  number of chromosomes: %s"%nchr)
    
    plt.figure(figsize = (15, 5))
    plt.title("Number of clusters estimation")
    plt.plot(np.arange(maxnumchr, 1, -1), mean_step_len[-maxnumchr+1:], 'b')
    plt.gca().invert_xaxis()
    plt.xlabel('number of clusters')
    plt.savefig(outbase+'.avg_step_len.svg', dpi=dpi)

    plt.figure()
    plt.title("Number of clusters estimation")
    plt.plot(np.arange(50, 1, -1), mean_step_len[-50+1:], 'b')
    plt.gca().invert_xaxis()
    plt.xlabel('number of clusters')
    plt.savefig(outbase+'.avg_step_len_50.svg', dpi=dpi)
                          
    logger(" Assigning contigs to chromosomes...")
    clusters = [[] for n in range(nchr)]
    assignments = sch.fcluster(Z, t=nchr, criterion='maxclust')
    for i, c in zip(assignments-1, bin_chr):
        clusters[i].append(c)
    return clusters
    
def bam2clusters(bam, fasta, outdir, minSize=2000, mapq=10, threads=4, dpi=100, upto=0, verbose=1, method="ward"):
    """Return clusters computed from from windowSizes"""
    logger("=== Clustering ===")
    outbase = os.path.join(outdir, "auto")
    
    # load clusters
    fname = outbase + ".clusters.tab_"
    if os.path.isfile(fname):
        clusters = [l[:-1].split('\t') for l in open(fname)]
        logger("  %s clusters loaded."%len(clusters))
        return clusters

    # get windows
    minSize, windows, chr2window, contig2size = contigs2windows(fasta, minSize, verbose)
    
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
        d, c2dists = bam2array(windows, contig2size, chr2window, bam,  mapq, upto, threads=threads, minSize=minSize)

        # save windows, array and plot
        logger("Saving array...")
        with gzip.open(outbase + ".windows.tab.gz", "w") as out:
            out.write("\n".join("\t".join(map(str, w)) for w in windows)+"\n")
        with open(outbase+".npz", "w") as out:
            np.savez_compressed(out, d)

        params = estimate_distance_parameters(outbase, contig2size=contig2size, c2dists=c2dists, windowSize=minSize)
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
        npy = np.load(outbase+".npz")
        d = npy[npy.files[0]]
        params = pickle.load(open(outbase+".distance.params"))
            
    # get clusters on transformed matrix
    #logger("Clustering...")
    #params = estimate_distance_parameters(outbase, bam, mapq, contig2size, minSize)
    transform = lambda x: distance_func(x+1, *params)
    clusters = cluster_contigs(outbase, transform(d), bin_chr, bin_position, dpi=dpi)
    
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
    parser.add_argument('--version', action='version', version='1.01e')   
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose")    
    parser.add_argument("-i", "--bam", nargs="+", help="BAM file(s)")
    parser.add_argument("-f", "--fasta", type=file, help="Genome FastA file")
    parser.add_argument("-o", "--outdir", required=1, help="output name")
    parser.add_argument("-m", "--minSize", default=2000, type=int,
                        help="minimum contig length (N90 or larger) [%(default)s]")
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
    