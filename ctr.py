#!/usr/bin/env python

import os, sys, gzip, subprocess, pysam
import numpy as np
from collections import Counter
from FastaIndex import FastaIndex

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

def bam2matches(bam, regions, mapq=10, upto=1e5):
    """Parse pysam alignments and return matching windows"""
    data = [0, 0]
    proc = _get_samtools_proc(bam, mapq, regions)
    # process reads from given library
    for i, line in enumerate(proc.stdout, 1):
        if upto and i>upto: break
        # unload data
        ref1, start1, mapq, cigar, ref2, start2, insertSize, seq = line.split('\t')[2:10]
        #start1, start2, seqlen = int(start1), int(start2), len(seq)
        # update ref if alignment in the same chrom
        if ref2 == "=":
            data[0] += 1
        else:
            data[1] += 1
    proc.terminate()
    return data

def main(bam):

    # get contig2size
    #faidx = FastaIndex(fasta)
    #contig2size = {c: faidx.id2stats[c][0] for c in faidx}
    sam = pysam.Samfile(bam)
    contig2size = {c: s for c, s in zip(sam.references, sam.lengths)}

    # estimate on largest contigs
    longest_contigs = sorted(contig2size, key=lambda x: contig2size[x], reverse=1)
    totdata = [0, 0]
    for c in longest_contigs[:25]:
        data = bam2matches(bam, regions=[c], mapq=10)
        print c, contig2size[c], 1.*data[0] / sum(data), sum(data)
        for i in range(2):
            totdata[i] += data[i]
    data = totdata
    print 1.*data[0] / sum(data), sum(data)

if __name__=='__main__':
    bam="/mnt/data/lpryszcz/cluster/hic/arath/idba/SRR2626163.contig.fa.bam"
    if len(sys.argv)>1:
        bam = sys.argv[1]
    main(bam)
        



