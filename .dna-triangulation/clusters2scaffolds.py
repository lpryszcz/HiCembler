#!/usr/bin/env python

import math, os, sys
import scipy.cluster.hierarchy as sch
import ete3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from collections import Counter
from datetime import datetime

import triangulation as tr
from array2clusters import array2clusters, get_name, get_shuffled

def normalize_rows(a):
    """Normalise rows of give array."""
    rows, cols = a.shape
    maxv = a.sum(axis=0).max()
    for i in xrange(rows):
        # only if any signal
        if a[i].max():
            a[i] *= maxv/a[i].sum()
    return a
    
def load_numpy(fname):
    """Load matrix from numpy dump"""
    # load array
    npy = np.load(fname)
    d = npy[npy.files[0]]
    # make symmetric
    d = d + d.T - np.diag(d.diagonal())
    d = normalize_rows(d)
    
    # load windows  
    windowfn = fname[:-4]+'.windows.tab.gz' 
    bin_chr = []
    bin_position = [] 
    for i, l in enumerate(gzip.open(windowfn)):
        chrom, start, end = l[:-1].split('\t')
        bin_chr.append(chrom)
        bin_position.append(map(int, (start, end)))
    return d, bin_chr, bin_position
    
def load_data_txt(file_txt, remove_nans=False, retain=1, chrs=None, remove_shorter=False, rename=False):

    '''
    load Hi-C interaction matrix from text file
    
    parameters:

    file_txt: file name. format "chr\tstart\tend\tdata1\tdata2\t..."
    remove_nans: removes nan rows/columns from all returned variables.
    retain: retain every x-th bin.
    chrs: load only these chrs. None mean load all chrs.
    
    returns:
    
    d: data matrix over the selected set of chromosomes.
    bin_chr: list of chr index assignment of each bin.
    bin_position: start and end position of each bin

    '''
    # NUMPY
    if file_txt.endswith(('.npz', '.npy')):
        d, bin_chr, bin_position = load_numpy(file_txt)
    # TEXT
    # open file, allowing gz
    elif file_txt.endswith('.gz'):
        fh = gzip.open(file_txt, 'r')
        d, bin_chr, bin_position = _load_raw(fh)
    else:
        fh = open(file_txt, 'r')
        d, bin_chr, bin_position = _load_raw(fh)

    # chromosome array
    bin_position = np.array(bin_position)
    bin_chr = np.array(bin_chr)

    # keep only relevant chromosomes
    if chrs:
        relevant_indices = np.any(bin_chr[None].T == chrs, 1)
        d = d[:, relevant_indices][relevant_indices, :]
        bin_chr = bin_chr[relevant_indices]
        bin_position = bin_position[relevant_indices, :]
        
    # rename chr to uniq names
    if rename:
        _bin_position = [(e-s) for s, e in bin_position]
        _bin_chr = ["%s_%s-%s" % (c, s, e) for c, (s, e) in zip(bin_chr, bin_position)]
        bin_position = np.array(bin_position)
        bin_chr = np.array(_bin_chr)
        
    # select subset of chromosomes / windows
    if retain > 1:
        d = d[::retain, ::retain]
        bin_chr = bin_chr[::retain]
        bin_position = bin_position[::retain, :]
        
    # eliminate nanas
    if remove_nans:
        valid_rowcols = ~(np.sum(np.isnan(d), 0) == d.shape[0])
        d = d[:, valid_rowcols][valid_rowcols, :]
        bin_chr = bin_chr[valid_rowcols]
        bin_position = bin_position[valid_rowcols, :]
        
    # eliminate zeros
    if remove_shorter:
        c = Counter(np.diff(bin_position, axis=1)[:, 0])
        windowSize, occurencies = c.most_common(1)[0]
        sys.stderr.write(" most common window: %s bp [%5.2f%s]\n"%(windowSize, occurencies*100./len(bin_chr), '%'))
        valid_rowcols = ~(np.diff(bin_position, axis=1)[:, 0]!=windowSize)
        d = d[:, valid_rowcols][valid_rowcols, :]
        bin_chr = bin_chr[valid_rowcols]
        bin_position = bin_position[valid_rowcols, :]   
    return d, bin_chr, bin_position


def clusters2scaffolds(infile, iterations=20, pnum=4,  evaluate=1, reduce_chr=1):
    """Compute scaffold for each cluster"""
    clustersFn = infile+".clusters.tab"
    if not os.path.isfile(clustersFn):
        tr.logger("Computing clusters...")
        clusters = array2clusters(infile)
    else:
        tr.logger("Loading precomputed clusters...")
        clusters = [l[:-1].split('\t') for l in open(clustersFn)]
    tr.logger(" loaded %s clusters."%len(clusters))
    
    # load matrix
    tr.logger("Loading matrix from %s ..."%infile)
    d, bin_chr, bin_position = tr.load_data_txt(infile, remove_nans=True, chrs=[], retain=1, remove_shorter=0)
    genomeSize = np.diff(bin_position, axis=1).sum()
    contig2size = {get_name(c): 0 for c in np.unique(bin_chr)}
    for c, (s, e) in zip(bin_chr, bin_position):
        contig2size[get_name(c)] += e-s
    print " loaded %s contigs summing %s bp"%(d.shape[0], genomeSize)

    #transform = lambda x: np.log(np.max(x+1))-np.log(x+1)
    #d = transform(d)
    
    # average contigs that share the same id
    if not evaluate and reduce_chr:
        logger("averaging contigs that share the same id...")
        d = tr.average_reduce_2d(d, bin_chr)
        np.unique(bin_chr)
              
    if evaluate:
        fig = plt.figure()
        mpl.rcParams['figure.subplot.hspace'] = 0.5
        mpl.rcParams['axes.titlesize'] = 10
        mpl.rcParams['axes.labelsize'] = 8
        mpl.rcParams['xtick.labelsize'] = 7
        mpl.rcParams['ytick.labelsize'] = 7
        x = y = int(math.sqrt(len(clusters)))
        if x*y < len(clusters):
            y += 1
            if x*y < len(clusters):
                x += 1
            
    tr.logger("Scaffolding %s clusters..."%len(clusters))
    for i, contigs in enumerate(clusters, 1):
        # get scaffold
        relevant_indices = np.any(bin_chr[None].T == contigs, 1)
        _d, _bin_position = d[:, relevant_indices][relevant_indices, :], bin_position[relevant_indices]      
        name = "cluster_%s"%i
        totsize = sum(e-s for s, e in _bin_position)
        sys.stderr.write(" %s %s %s kb in %s contigs\n"%(i, name, totsize/1000, _d.shape[0]))
        scales, pos, x0, fvals = tr.assemble_chromosome(_d, pnum=pnum, iterations=iterations, shuffle=True, return_all=True)
        # how to correlate estimated position with real position?
        # plot
        if evaluate:
            ax = fig.add_subplot(x, y, i)
            ax.set_title(name)
            plt.plot(_bin_position, pos[0, :], 'b.')
            # plot axes labels only on edges
            if i >= len(clusters)-x:
                plt.xlabel("Expected position")
            if i%y==1:
                plt.ylabel("Predicted position")
    if evaluate:
        tr.logger("Saving figure...")
        fig.savefig(infile+'.pred_position.svg')
    tr.logger("Done!")
    
def main():
    import argparse
    parser  =  argparse.ArgumentParser(description='De novo assembly of Hi-C data.', 
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--infile', required=True,  
                        help='Hi-C interaction matrix input file')
    parser.add_argument('-iter', '--iterations', default=10, type=int, 
                        help='no. of iterations [%(default)s]')
    parser.add_argument('-t', '--threads', default=4, type=int, 
                        help='no. of threads to use [%(default)s]')
    parser.add_argument('-e', '--evaluate', action='store_true', default=False, 
                        help='evaluation mode. chromosome names are assumed to be the true chromosomal assignment.')
    
    o = parser.parse_args()
         
    #infile = sys.argv[1] 
    clusters2scaffolds(o.infile, o.iterations, o.threads, o.evaluate)
    
if __name__=="__main__":
    t0 = datetime.now()
    main()
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
