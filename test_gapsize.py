#!/usr/bin/env python

#from array2scaffolds import load_matrix
import os, sys, gzip
import numpy as np
import matplotlib.pyplot as plt#; plt.ion()
from scipy.optimize import curve_fit
from datetime import datetime
from collections import Counter
from bam2scaffolds import bam2array, normalize_rows
import pysam

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

def get_name(contig):
    return contig.split()[0] 

def normalize_rows(a):
    """Normalise rows so the sums among rows are identical."""
    rows, cols = a.shape
    maxv = a.sum(axis=0).max()
    for i in xrange(rows):
        # only if any signal
        if a[i].max():
            a[i] *= 1.*maxv/a[i].sum()
    return a

def get_contig2size(bin_chr, bin_position):
    """Return contig2size"""
    # calculate genome size and contig2size
    contig2size = {get_name(c): 0 for c in np.unique(bin_chr)}
    for c, (s, e) in zip(bin_chr, bin_position):
        contig2size[get_name(c)] += e-s
    return contig2size
    
def load_matrix(fname, chrs=[], remove_shorter=True, scaffolds=[], verbose=0):
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
        
    # make symmetric at the end to save memory and time
    d += d.T
    d -= np.diag(d.diagonal()/2)
    d = normalize_rows(d)
    d = normalize(d)
    #d, bin_chr, bin_position = normalize_diagional(d, bin_chr, bin_position)
    return d, bin_chr, bin_position, contig2size

t0 = datetime.now()
d2c = {}

'''
print "Loading..."
d, bin_chr, bin_position, contig2size = load_matrix('/home/lpryszcz/cluster/hic/arath/_archives/bam2scaffolds_new/10k.npz')

# estimate on largest contigs
longest_contigs = sorted(contig2size, key=lambda x: contig2size[x], reverse=1)

dists, contacts = [], []
# get middle of each window
positions = bin_position[:,0] + (np.diff(bin_position, axis=1)/2)[:,0]

print "Populating..."
for c in longest_contigs[:1]:
  indices = np.argwhere(bin_chr==c)[:,0][:1000]
  # add diagonal
  d2c[0] = [d[i][i] for i in indices]
  for ii, i in enumerate(indices[:-1]): #d[np.triu_indices(indices.shape[0], 1)]
    for j in indices[ii+1:]:
      dist = (positions[j]-positions[i]) / 1000
      if dist not in d2c:
        d2c[dist] = []
      d2c[dist].append(d[i][j])
      #dists.append(dist)
      #contacts.append(d[i][j])
'''

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

def main(bam, windowSize=2000, mapq=0):
    d2c[0] = []

    print "Loading..."
    sam = pysam.Samfile(bam)
    contig2size = {c: s for c, s in zip(sam.references, sam.lengths)}
    # estimate on largest contigs
    longest_contigs = sorted(contig2size, key=lambda x: contig2size[x], reverse=1)
    #c = longest_contigs[0]
    for c in longest_contigs[:5]:
      s = contig2size[c]
      sys.stderr.write(' %s %s bp   \r'%(c, s))
      n = s / windowSize + 1
      positions = range(windowSize, n*windowSize, windowSize)
      arrays = [np.zeros((n,n), dtype='float32')]
      chr2window = [{c: 0}]
      arrays = bam2array(arrays, [windowSize], chr2window, [bam], mapq, regions=[c], upto=1e6, verbose=0)
      d = arrays[0]#; print c, len(positions), s, d.shape
      
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
    #'''

    dists = np.array(sorted(d2c)[:50])
    contacts = [d2c[d] for d in dists]

    print "Plotting %s distances..."%len(contacts)
    plt.title("HiC contacts at given distance")
    plt.boxplot(contacts, 1, '', positions=dists, widths=.75*dists[1])#; plt.legend("HiC data")
    plt.xticks(rotation=90)

    plt.xlabel("distance [kb]")
    plt.ylabel("contacts")
    plt.xlim(xmin=-dists[1])
    plt.ylim(ymin=0)#, ymax=20)
    #plt.yscale('log')

    print "fitting curve..." # http://stackoverflow.com/a/11209147/632242
    x = dists; yn = np.array([np.median(c) for c in contacts])#, dtype='double') [2:]

    step = dists[1]/4.
    xs = np.arange(0, max(x)+step, step)
    # Non-linear fit
    #def func(x, b, c, d): return np.exp(-b * x) / c + 1 / (d*x)
    def func(x, d, e): return 1 / (d*x**e)
    #def func(x, d): return 1 / (d*x)
    popt, pcov = curve_fit(func, x[5:], yn[5:]); print pcov
    plt.plot(xs[1:], func(xs[1:], *popt), 'r-', label="Non-linear fit\n$ 1 / {%0.2f x^ {%0.2f}} $"%tuple(popt))

    def ddd(x, b, c): return np.exp(-b * x) / c 
    #def ddd(x, b, c, d): return np.exp(-b * x) / c + d
    # def ddd(x, b, c, d): return np.exp(-b * x) / c + 1 / d
    popt2, pcov2 = curve_fit(ddd, x, yn); print pcov2
    plt.plot(xs, ddd(xs, *popt2), 'b--', label="Distance dependent decay\n$ y = e^{-%0.2f x} / %0.5f $"%tuple(popt2))

    plt.legend(fancybox=True, shadow=True)
    plt.savefig(bam+".fit.png")
    plt.ylim(ymin=0, ymax=yn[0]/30.)
    plt.savefig(bam+".fit.zoom.png")
    dt = datetime.now() - t0; print dt
    
if __name__=='__main__':
    bam = "/home/lpryszcz/cluster/hic/arath/platanus/ARATH.d05.l100_contig.fa.bam"
    windowSize = 2000
    if len(sys.argv)>1: 
        bam = sys.argv[1]
    if len(sys.argv)>2: 
        windowSize = int(sys.argv[2])
    main(bam, windowSize)


