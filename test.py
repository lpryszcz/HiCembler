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

def get_longest(t, maxdist=6):
    """Return node having longest branch"""
    #n = sorted(t.traverse(), key=lambda n: 2*n.dist-t.get_distance(n), reverse=1)[0]
    n = t
    bestdist = 2*n.dist-n.get_distance(t)
    for _n in t.traverse():
        if _n.get_distance(t, topology_only=1) > maxdist:
            break
        if 2*_n.dist-_n.get_distance(t) > bestdist:
            n = _n
            bestdist = 2*_n.dist-_n.get_distance(t)
    return n
    
def get_subtrees(d, bin_chr, bin_position, method="ward", nchrom=1000, distfrac=0.75):
    maxtdist = 0
    i = 0
    subtrees = []
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
          
        # get longest branch
        n = get_longest(t)
        pruned = n.get_leaf_names()         
        subtrees.append(pruned)
        # prune array        
        indices = [_i for _i, name in enumerate(names) if name not in set(pruned)]
        d = d[indices, :]
        d = d[:, indices]
        bin_chr = bin_chr[indices]
        bin_position = bin_position[indices, :]
            
    if i:    
        subtrees.append(t.get_leaf_names())
    return subtrees
    

def main(fn):
    d, bin_chr, bin_position, contig2size = load_matrix(fn, remove_shorter=0)
    d = transform(d)
    subtrees = get_subtrees(d, bin_chr, bin_position)
    
    logger(" Assigning contigs to %s clusters..."%len(subtrees))
    total = correct = 0
    contig2cluster = {get_name(c): Counter() for c in np.unique(bin_chr)}
    for i, subtree in enumerate(subtrees, 1):
        c = Counter(map(get_name, subtree))
        total += len(subtree)
        correct += c.most_common(1)[0][1]
        # poplate contig2clustre
        for k, v in c.iteritems():
            if not k: continue
            contig2cluster[get_name(k)][i] += v
    logger("  %s / %s [%.2f%s]"%(correct, total, 100.*correct/total, '%'))

    logger(" Weak assignments...")
    clusters = [[] for i in range(len(subtree))]
    withoutCluster, weakCluster = [], []
    for c, counter in contig2cluster.iteritems():
        if not counter:
            withoutCluster.append(c)
            continue
        # get major cluster
        clusteri, count = counter.most_common(1)[0]
        mfrac = 1. * count / sum(counter.itervalues())
        clusters[clusteri].append(c)
        if mfrac<.66:
            weakCluster.append(c)
    logger("  %s bp in %s contigs without assignment."%(sum(contig2size[c] for c in withoutCluster), len(withoutCluster)))
    logger("  %s bp in %s contigs having weak assignment."%(sum(contig2size[c] for c in weakCluster), len(weakCluster)))
      
    outfile = fn[:-4]+".clusters.tab"
    clusters = filter(lambda x: x, clusters)
    totsize = 0
    logger("Reporting %s clusters to %s ..."%(len(clusters), outfile))
    with open(outfile, "w") as out:
        for i, cluster in enumerate(clusters, 1):
            #print " cluster_%s %s windows; %s"%(i, len(cluster), Counter(get_chromosome(cluster).most_common(3)))
            clSize = sum(contig2size[c] for c in cluster)
            totsize += clSize
            out.write("\t".join(cluster)+"\n")
    logger("  %3s bp in %s clusters."%(totsize, len(clusters)))
        
if __name__=="__main__":
    t0 = datetime.now()
    fn = '/home/lpryszcz/cluster/hic/arath/_archives/snap/SRR2626163.100k.npz'
    if len(sys.argv)>1:
      fn = sys.argv[1]
    main(fn)
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)


