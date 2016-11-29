#!/usr/bin/env python

import os, sys, gzip
import numpy as np
from collections import Counter
import fastq2array as fa

def get_name(contig):
    return contig.split()[0].split('|')[-1].split('.')[0]

def get_contig2size(bin_chr, bin_position):
    """Return contig2size"""
    # calculate genome size and contig2size
    contig2size = {get_name(c): 0 for c in np.unique(bin_chr)}
    for c, (s, e) in zip(bin_chr, bin_position):
        contig2size[get_name(c)] += e-s
    return contig2size

def main(fname):
    # load array
    npy = np.load(fname)
    d = npy[npy.files[0]]
    d += d.T - np.diag(d.diagonal())
    print d.max(), d.mean()
    #d = d.T
    #d /= d.mean(); print d.max(), d.mean()

    verbose, dpi, ext = 1, 300, "png"

    windowfn = fname[:-4]+'.windows.tab.gz'
    bin_chr = []
    bin_position = [] 
    for i, l in enumerate(gzip.open(windowfn)):
        chrom, start, end = l[:-1].split('\t')
        bin_chr.append(get_name(chrom))
        bin_position.append(map(int, (start, end)))
        
    # chromosome array
    bin_position = np.array(bin_position)
    bin_chr = np.array(bin_chr)
    contig2size = get_contig2size(bin_chr, bin_position)

    c = Counter(np.diff(bin_position, axis=1)[:, 0])
    windowSize, occurencies = c.most_common(1)[0]; print windowSize, occurencies, len(bin_chr)

    base2chr = {}
    genomeSize = 0
    for c in np.unique(bin_chr):
        base2chr[genomeSize+contig2size[c]/2] = c
        genomeSize += contig2size[c]

    print "Saving figure as %s.%s" % (fname, ext)
    fa.plot(fname, d, genomeSize, base2chr, windowSize, dpi, ext=ext)

if __name__=='__main__':
    fnames = ["/home/lpryszcz/cluster/hic/arath/_archives/bam2scaffolds/100k.npz"]
    if len(sys.argv)>1:
        fnames = sys.argv[1:]; print fnames
    for fname in fnames:
        main(fname)


