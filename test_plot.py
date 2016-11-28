#!/usr/bin/env python

import numpy as np
import scipy.cluster.hierarchy as sch
from sklearn.cluster import ward_tree, AffinityPropagation, MeanShift, DBSCAN, Birch, KMeans
import ete3, gzip, os, resource, sys
#from array2scaffolds import load_matrix, logger, transform
from collections import Counter
from datetime import datetime
import fastcluster
import matplotlib.pyplot as plt

def logger(message, log=sys.stdout):
    """Log messages"""
    memory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    log.write("[%s] %s    [memory: %6i Mb]\n"%(datetime.ctime(datetime.now()), message, memory))

transform = lambda x: np.log(np.max(x+1))-np.log(x+1)

# update sys.path & environmental PATH
root = os.path.dirname(os.path.abspath(sys.argv[0]))
src = ["bin", "bin/snap", "bin/sinkhorn_knopp"]
paths = [os.path.join(root, p) for p in src]
sys.path = paths + sys.path
os.environ["PATH"] = "%s:%s"%(':'.join(paths), os.environ["PATH"])

from sinkhorn_knopp import sinkhorn_knopp

def normalize(d, bin_chr, bin_position, max_iter=1000, epsilon=0.0001, windowSize=1000.):
    """Return symmetric and fully balanced matrix using SinkhornKnopp"""
    # make symmetric & normalise
    d += d.T
    d -= np.diag(d.diagonal()/2)
    #return d, bin_chr, bin_position
    # normalize by windows size
    sizes = np.diff(bin_position, axis=1)#[:, 0]
    #c = Counter(sizes.reshape(len(sizes)))
    #windowSize, occurencies = c.most_common(1)[0]; print windowSize, occurencies
    #d *= 1. * windowSize / sizes
    #d *= windowSize **2 / (sizes*sizes.T)**0.5; print sizes.shape, sizes.T.shape #reshape(len(sizes),1))
    # full balancing
    sk = sinkhorn_knopp.SinkhornKnopp(max_iter=max_iter, epsilon=epsilon); d += 1; d /= d.max(); d = sk.fit(d) #* 100000
    # 1 round balancing
    #sk = sinkhorn_knopp.SinkhornKnopp(max_iter=1); d += 1; d /= d.max(); d = sk.fit(d)
    '''
    axis = 1; d *= 1. * d.sum(axis=axis).max() / d.sum(axis=axis); print "axis %s norm"%axis #normalize_rows(d)
    ''' # diagonal mean normalisation
    # normalize_rows(d)
    '''
    indices = d.diagonal()!=0; print "diag norm"
    d = d[indices, :]
    d = d[:, indices]
    bin_chr = bin_chr[indices]
    bin_position = bin_position[indices, :]    
    d *= np.mean(d.diagonal()) / d.diagonal() #'''
    return d, bin_chr, bin_position

def get_contig2size(bin_chr, bin_position):
    """Return contig2size"""
    # calculate genome size and contig2size
    contig2size = {get_name(c): 0 for c in np.unique(bin_chr)}
    for c, (s, e) in zip(bin_chr, bin_position):
        contig2size[get_name(c)] += e-s
    return contig2size

def load_matrix(fname, chrs=[], remove_shorter=True, scaffolds=[], verbose=0, remove_nans=1, remove_zeros=1):
    """Load Hi-C interaction matrix from numpy dump
    generated by fastq2array.py. 
     
    Returns:
    d: data matrix over the selected set of chromosomes.
    bin_chr: list of chr index assignment of each bin.
    bin_position: start and end position of each bin
    """
    if scaffolds:
        remove_shorter = True
        
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
    contig2size = get_contig2size(bin_chr, bin_position)
    
    # eliminate nanas
    if remove_nans:
        indices = ~(np.sum(np.isnan(d), 0) == d.shape[0])
        if indices.sum() < d.shape[0]:
            print "remove_nans:", indices.sum(), d.shape
            d = d[indices, :]
            d = d[:, indices]
            bin_chr = bin_chr[indices]
            bin_position = bin_position[indices, :]
            
    if remove_zeros:
        indices = ~(np.any((np.sum(d, axis=0)==0, np.sum(d, axis=1)==0), axis=0))
        if indices.sum() < d.shape[0]:
            print "removed rows/columns summing to zero:", indices.sum(), d.shape
            d = d[indices, :]
            d = d[:, indices]
            bin_chr = bin_chr[indices]
            bin_position = bin_position[indices, :]
    
    #''' # eliminate 
    c = Counter(np.diff(bin_position, axis=1)[:, 0])
    windowSize, occurencies = c.most_common(1)[0]    
    if remove_shorter:
        if verbose:
            sys.stderr.write(" most common window: %s bp [%5.2f%s]\n"%(windowSize, occurencies*100./len(bin_chr), '%'))
        indices = ~(np.diff(bin_position, axis=1)[:, 0]!=windowSize)
        d = d[indices, :]
        d = d[:, indices]
        bin_chr = bin_chr[indices]
        bin_position = bin_position[indices, :]
    #'''
    else:
        # normalise by length
        sizenorm = np.array([1.0*windowSize/(e-s) for s, e in bin_position])
        d *= sizenorm #'''
   
    # keep only relevant chromosomes
    if chrs:
        indices = np.any(bin_chr[None].T == chrs, 1)
        d = d[indices, :]
        d = d[:, indices]
        bin_chr = bin_chr[indices]
        bin_position = bin_position[indices, :]

    # combine existing array using information from previous round of scaffolding
    if scaffolds: 
        contig2indices = get_contig2indices(bin_chr)
        indices, bin_chr, bin_position = [], [], []
        for i, scaffold in enumerate(scaffolds, 1):
            name = "scaffold%s"%i
            indices += get_indices(scaffold, contig2indices)
            bin_chr += [name]*len(indices)
            bin_position += [(s, s+windowSize) for s in range(0, windowSize*len(indices), windowSize)]
        # combine
        d = d[:, indices][indices, :]
        bin_chr = np.array(bin_chr)
        bin_position = np.array(bin_position)
        contig2size = get_contig2size(bin_chr, bin_position)
    
    d, bin_chr, bin_position = normalize(d, bin_chr, bin_position)
            
    return d, bin_chr, bin_position, contig2size

def get_names(bin_chr, bin_position):
    return ["%s %s"%(c, s) for c, (s, e) in zip(bin_chr, bin_position)]

def get_name(contig):
    return contig.split()[0]
    
def get_chr_name(n):
    return n.split()[0].split(".")[0]
    
def get_chromosome(names): return Counter(get_chr_name(n) for n in names)    


def main(fn, method="ward"): #
    d, bin_chr, bin_position, contig2size = load_matrix(fn, remove_shorter=0)
    sizes = np.diff(bin_position, axis=1)[:, 0] / 1000
    contacts = d.diagonal()
    print d.sum(), d.diagonal().sum()
    # get bins
    bins = np.arange(1, 101, 5)
        
    # get counts
    contacts = [[] for i in range(len(bins)+1)]
    for s, c in zip(np.digitize(sizes, bins, right=1), d.diagonal()):
        contacts[s].append(c)
    print len(contacts), len(bins), len(sizes)#, contacts# np.digitize(sizes, bins, right=1)
    plt.title("HiC contacts at given distance")
    plt.boxplot(contacts[:-1], 1, '', positions=bins, widths=.75*bins[0])#; plt.legend("HiC data")
    plt.xticks(rotation=90)

    plt.xlabel("contig size [kb]")
    plt.ylabel("self contacts")
    plt.xlim(xmin=-bins[0])
    plt.ylim(ymin=0)#, ymax=20)
    #plt.yscale('log')
    #plt.show()
    outfn = fn+".selfcontacts.png"
    plt.savefig(outfn)
    print "Figure saved as: %s"%outfn
        
if __name__=="__main__":
    t0 = datetime.now()
    method = "ward"
    fn = '/home/lpryszcz/cluster/hic/arath/_archives/snap/SRR2626163.100k.npz'
    if len(sys.argv)>1:
      fn = sys.argv[1]
    if len(sys.argv)>2:
      method = sys.argv[2]
    main(fn, method)
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)

