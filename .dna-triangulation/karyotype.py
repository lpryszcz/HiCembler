#!/usr/bin/env python

import argparse
import scipy as sp
import scipy.cluster
import sys
import collections
import triangulation as tr
import numpy as np
import matplotlib.pyplot as plt

from triangulation import logger
from datetime import datetime

def karyotype(infile, outfile, nchr, drop, included_chrs, seed, rand_frac, rand_n, evaluate, maxnumchr=1000): 
    """Estimate chromosome number & evaluate"""
    logger("loading matrix...")
    d, bin_chr, bin_position = tr.load_data_txt(infile, remove_nans=True, chrs=[], retain=drop, remove_shorter=0)
    ncontigs = bin_chr.shape[0]
    genomeSize = np.diff(bin_position, axis=1).sum()
    logger(" loaded %s contigs summing %s bp"%(ncontigs, genomeSize))
    # adjust maxnumchr to avoid errors
    if ncontigs < maxnumchr*2:
        maxnumchr = ncontigs/2
        sys.stderr.write("  adjusted maxnumchr to %s\n"%maxnumchr)

    # get chromosome names if not provided
    if not included_chrs:
        starts = ("ENA|", "gi|","gb|")
        chrnames = lambda x: x.startswith('chr') and len(x)<10 or x.startswith(starts)
        included_chrs = filter(chrnames, set(bin_chr)) # chrXIII

    logger("karyotyping...")
    pred_nchr = False
    if nchr == 0:
        nchr = maxnumchr
        pred_nchr = True

    n = d.shape[0]
    transform = lambda x: np.log(np.max(x+1))-np.log(x+1)
    res = tr.predict_karyotype(d, nchr=nchr, pred_nchr=pred_nchr, transform=transform, shuffle=0, #True, 
                               seed=seed, rand_frac=rand_frac, rand_n=rand_n)
    if pred_nchr:
        clust, Z, nchr, mean_step_len, wrong = res
        if wrong:
            bin_chr = np.delete(bin_chr, wrong, 0)
            bin_position= np.delete(bin_position, wrong, 0)
            n -= len(wrong)
        logger(" identified %s chromosomes."%nchr)
        
        np.savetxt(outfile+'_avg_step_len.tab', np.c_[np.arange(maxnumchr, 1, -1), mean_step_len[-maxnumchr+1:]], fmt='%s', delimiter='\t')
        np.savetxt(outfile+'_clusteringZ.tab', Z, fmt='%s', delimiter='\t')
        np.savetxt(outfile+'_clusters.tab', np.c_[bin_chr, bin_position, clust], fmt='%s', delimiter='\t')
        
        logger(" plotting...")
        plt.figure(figsize = (15, 5))
        plt.plot(np.arange(maxnumchr, 1, -1), mean_step_len[-maxnumchr+1:], 'b')
        plt.gca().invert_xaxis()
        plt.xlabel('number of clusters')
        plt.savefig(outfile+'_avg_step_len.svg', dpi=600)

        plt.figure()
        plt.plot(np.arange(80, 1, -1), mean_step_len[-80+1:], 'b')
        plt.gca().invert_xaxis()
        plt.xlabel('number of clusters')
        plt.savefig(outfile+'_avg_step_len_80.svg', dpi=600)

        sys.setrecursionlimit(100000)
        #tr.plot_dendro(outfile+"_dendro.svg", Z)
    else:
        clust, Z, wrong = res
        if wrong:
            bin_chr = np.delete(bin_chr, wrong, 0)
            bin_position= np.delete(bin_position, wrong, 0)
            n -= len(wrong)

    if evaluate and included_chrs:
        logger("evaluating...")
        # match each cluster to the chromosome which most of its members belongs to
        chr_order = dict(zip(included_chrs, range(len(included_chrs))))
        new_clust = np.zeros(n, dtype=bin_chr.dtype)
        new_clust_num = np.nan*np.ones(n)
        for i in range(nchr):
            chrname = collections.Counter(bin_chr[clust == i]).most_common(1)[0][0]
            # make sure all chromosomes are present in reference
            if chrname in chr_order:
                new_clust[clust == i] = collections.Counter(bin_chr[clust == i]).most_common(1)[0][0]
                new_clust_num[clust == i] = chr_order[chrname]
            
        # calculate accuracy
        accuracy = np.sum(new_clust == bin_chr)/float(n)
        logger(" estimated accuracy: %.5f"%accuracy)
        # plot figure
        plt.figure(figsize = (15, 5))
        tr.chr_color_plot(np.mean(bin_position, 1), bin_chr, new_clust_num, included_chrs, int(genomeSize*0.001))
        plt.savefig(outfile+'_evaluation.svg', dpi=600)
        np.savetxt(outfile+'_evaluation.tab', np.c_[bin_chr, bin_position, new_clust], fmt='%s', delimiter='\t')

def main():
    parser  =  argparse.ArgumentParser(description='De novo karyotyping of Hi-C data.', 
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-in', dest='infile', type=str, required=True,  
                        help='Hi-C interaction matrix input file')
    parser.add_argument('-out', type=str, required=True, dest='outfile',
                        help='prefix for output files')
    parser.add_argument('-nchr', dest='nchr', type=int, default=0, 
                        help='number of chromosomes/clusters. 0 will automatically estimate this number.')
    parser.add_argument('-drop', dest='drop', type=int, default=1, 
                        help='leaves every nth bin in the data,  ignoring the rest. 1 will use whole dataset.')
    parser.add_argument('-ci', dest='included_chrs', nargs='*', type=str, default=[], 
                        help='list of chromosomes/contigs to include. If empty,  uses all chromosomes.')
    parser.add_argument('-s', dest='seed', type=int, default=0, 
                        help='seed for randomizations')
    parser.add_argument('-f', dest='rand_frac', type=float, default=0.8, 
                        help='fraction of data to use for average step length calculation')
    parser.add_argument('-n', dest='rand_n', type=int, default=20, 
                        help='number of iterations for average step length calculation')
    parser.add_argument('-e', dest='evaluate', action='store_true', 
                        help='evaluation mode. chromosome names are assumed to be the true chromosomal assignment.')
    
    o = parser.parse_args()
    
    karyotype(o.infile, o.outfile, o.nchr, o.drop, o.included_chrs, o.seed,
              o.rand_frac, o.rand_n, o.evaluate)
    
if __name__ == "__main__":
    t0 = datetime.now()
    main() 
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
