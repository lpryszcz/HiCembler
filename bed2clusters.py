#!/usr/bin/env python

import os, sys
import numpy as np
from datetime import datetime

def get_name(contig):
    return contig.split()[0].split('|')[-1].split('|')[0]

def test():
    out = sys.stdout
    fasta, ref = sys.argv[1:3]
    
    # generate & load contig2chrom
    if not os.path.isfile("%s.bed"%fasta):
        # generate index
        if not os.path.isfile("%s.suf"%ref):
            os.system("lastdb %s %s"%(ref, ref))
        # generate chromosome to tab
        os.system("lastal -l 100 -C 2 -P 4 %s %s | last-split - | maf-convert tab - | tab2chromosome.py > %s.bed"%(ref, fasta, fasta))
        
    c2chr = {l.split('\t')[3]: get_name(l.split('\t')[0]) for l in open("%s.bed"%fasta)}
    
    chr2contigs = {}
    for c, chrom in c2chr.iteritems():
      if chrom not in chr2contigs:
        chr2contigs[chrom] = set()
      chr2contigs[chrom].add(c)
    
    clusters = chr2contigs.values()
    for cluster in clusters:
      out.write('\t'.join(cluster)+'\n')
    
            
if __name__=="__main__":
    t0 = datetime.now()
    #main()
    test()
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)


