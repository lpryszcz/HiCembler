#!/usr/bin/env python

import numpy as np; import os, resource
from array2scaffolds import *
from collections import Counter
from datetime import datetime

def get_chromosome(names): return Counter(n.split(" ")[0].split(".")[0] for n in names)

def truncate(t, mind=3, maxd=0):
    for i, n in enumerate(t.traverse(), 1):
        dist = t.get_distance(n, topology_only=1)
        chrs = get_chromosome(n.get_leaf_names())        
        if dist>mind and len(chrs)==1 or maxd and dist>maxd:
            n.leaves = n.get_leaf_names()
            n.chrs = chrs
            n.name="%s %s chrs %s leaves"%(chrs.most_common(1)[0][0], len(chrs), len(n))
            for _n in n.get_children():
                n.remove_child(_n)
    return t

def get_names(bin_chr, bin_position):
    return ["%s %s"%(c, s) for c, (s, e) in zip(bin_chr, bin_position)]

def main2(fn, method="ward", nchrom=1000, distfrac=0.75): #ward average

    maxtdist = 0
    d, bin_chr, bin_position, contig2size = load_matrix(fn)
    d = transform(d)
    i = 0
    clusters = []
    for i in range(1, nchrom):
    
        Z = sch.linkage(d[np.triu_indices(d.shape[0], 1)], method=method)
        
        tree = sch.to_tree(Z, False)
        names = get_names(bin_chr, bin_position)
        nw = getNewick(tree, "", tree.dist, names)    
        t = ete3.Tree(nw)
    
        tname, tdist = t.get_farthest_leaf()#[1]
        if maxtdist < tdist:
            maxtdist = t.get_farthest_leaf()[1]
        # break if small subtree
        if tdist < maxtdist*distfrac:
            break
          
        print i, len(t), tname.name, tdist, get_chromosome(t.get_leaf_names())        
        # get longest branch
        dists = sorted((n.dist for n in t.traverse()), reverse=1)
        # 
        n = sorted(t.traverse(), key=lambda n: 2*n.dist-t.get_distance(n), reverse=1)[0]
        print " %s %s pruned; %s; %s"%(len(n), n.dist, get_chromosome(n.get_leaf_names()), dists[:10])

        pruned = n.get_leaf_names()         
        #'''
        # prune array        
        indices = [_i for _i, name in enumerate(names) if name not in set(pruned)]
        d = d[indices, :]
        d = d[:, indices]
        bin_chr = bin_chr[indices]
        bin_position = bin_position[indices, :]

        t2 = truncate(ete3.Tree(t.write())); t2.render("%s.pdf"%i)
        clusters.append(pruned)
            
    if i:    
      # remove last subcluster from tree
      clusters.append(t.get_leaf_names())
      #t.get_ancestors()[0].remove_child(n); 
      t = truncate(t); t.render("%s.pdf"%i)
    
    print "Reporting %s clusters..."%len(clusters)
    for i, cluster in enumerate(clusters, 1):
        print " cluster_%s %s windows; %s"%(i, len(cluster), Counter(get_chromosome(cluster)))


def main(fn, method="ward", nchrom=1000, distfrac=0.75): #ward average

    maxtdist = 80
    d, bin_chr, bin_position, contig2size = load_matrix(fn)
    d = transform(d)

    Z = sch.linkage(d[np.triu_indices(d.shape[0], 1)], method=method)
    
    tree = sch.to_tree(Z, False)
    names = get_names(bin_chr, bin_position)
    nw = getNewick(tree, "", tree.dist, names)    
    t = ete3.Tree(nw)
    maxtdist = t.get_farthest_leaf()[1]

    i = 0
    clusters = []
    for i in range(1, nchrom):
        tname, tdist = t.get_farthest_leaf()#[1]
        # break if small subtree
        if tdist < maxtdist*distfrac:
            break
          
        print i, len(t), tname.name, tdist, get_chromosome(t.get_leaf_names())        
        # get longest branch
        dists = sorted((n.dist for n in t.traverse()), reverse=1)
        # 
        n = sorted(t.traverse(), key=lambda n: 2*n.dist-t.get_distance(n), reverse=1)[0]
        print " %s %s pruned; %s; %s"%(len(n), n.dist, get_chromosome(n.get_leaf_names()), dists[:10])

        pruned = n.get_leaf_names()         

        t2 = truncate(ete3.Tree(t.write())); t2.render("%s.pdf"%i)
        # pruning is slower than recalculating the tree!! 18vs11s or 44 vs 69s
        t.prune([_n for _n in t.get_leaf_names() if _n not in set(pruned)], preserve_branch_length=1)
        #n.get_ancestors()[0].remove_child(n); n.delete()
        clusters.append(pruned)
            
    if i:    
      # remove last subcluster from tree
      clusters.append(t.get_leaf_names())
      #t.get_ancestors()[0].remove_child(n); 
      t = truncate(t); t.render("%s.pdf"%i)
    
    print "Reporting %s clusters..."%len(clusters)
    for i, cluster in enumerate(clusters, 1):
        print " cluster_%s %s windows; %s"%(i, len(cluster), Counter(get_chromosome(cluster)))

if __name__=="__main__":
    t0 = datetime.now()
    fn = '/home/lpryszcz/cluster/hic/arath/_archives/snap/SRR2626163.100k.npz'
    if len(sys.argv)>1:
      fn = sys.argv[1]
    main2(fn)
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)


