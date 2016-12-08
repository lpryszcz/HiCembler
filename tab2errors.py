#!/usr/bin/env python
# Align contigs onto chromosomes, parse the alignments (.tab.gz) and report synteny blocks. 
# The difference in gap size will be calculated, but only within each synteny block. 

import os, sys, gzip
from datetime import datetime
#from FastaIndex import FastaIndex
from collections import Counter
	
def generate_tab(ref, fasta, threads=4):
	"""Generate .tab.gz its fname"""
	if not os.path.isfile("%s.tab.gz"%fasta):
        # generate index
		if not os.path.isfile("%s.suf"%ref):
			os.system("lastdb %s %s"%(ref, ref))
        # generate chromosome to tab
		cmd0 = "lastal -l 100 -C 2 -P %s %s %s | last-split - | maf-convert tab - | gzip > %s.tab.gz"
		os.system(cmd0%(threads, ref, fasta, fasta))
	return "%s.tab.gz"%fasta

def get_best_chr(matches):
	"""Return best chr for given query matches"""
	total_alg = 0
	talg = Counter()
	for qstart, qend, t, tstart, tend, strand in matches:
		total_alg += tend-tstart
		talg[t] += tend-tstart
	bestt, talg = talg.most_common(1)[0]
	
	strandc = Counter()
	for qstart, qend, t, tstart, tend, strand in matches:
		if t!=bestt: 
			continue
		strandc[strand] += tend-tstart
	beststrand, salg = strandc.most_common(1)[0]#; print strandc
	return bestt, talg, beststrand, total_alg

def get_blocks(matches, bestt, beststrand):
	"""Parse matches and report synteny blocks"""
	qstart, qend, t, tstart, tend, strand = matches[0]
	blocks = [[(t, tstart, tend-tstart, strand),]]#matches[0][2:6],]]
	event = False
	reverse = False
	for i in xrange(1, len(matches)):
		qstart, qend, t, tstart, tend, strand = matches[i]#; print matches[i]
		# inversion or inter-chromosomal translocation 
		# or intra-chromosomal translocation: tstart < previos for + or tstart > previous for -
		if strand != matches[i-1][-1] or t != matches[i-1][2] \
		  or strand == "+" and tstart < matches[i-1][3]\
		  or strand == "-" and tstart > matches[i-1][3]:
			blocks.append([])
		elif blocks[-1]:
			# get gap size difference
			estgap = abs(qstart-matches[-1][1])
			trugap = abs(tstart-matches[-1][4])
		# add match
		blocks[-1].append((t, tstart, tend-tstart, strand))#matches[i][2:6])
	return blocks

def get_name(name): return name.split()[0].split('|')[-1].split('.')[0]
	
def process_hits(ref, fasta, identityTh=0.9, overlapTh=.66, minscore=100, threads=4):
    """Generate & process LASTal hits"""
    # execute last and get best query-to-reference matches only
    tabfn = generate_tab(ref, fasta, threads)
    q2matches = {}
    queries = {}
    targets = {}
    for l in gzip.open(tabfn): 
        if l.startswith('#'):
            continue
        # unpack
        (score, t, tstart, talg, tstrand, tsize, q, qstart, qalg, qstrand, qsize, blocks) = l.split()[:12]
        (score, qstart, qalg, qsize, tstart, talg, tsize) = map(int, (score, qstart, qalg, qsize, tstart, talg, tsize))
        if score<minscore:
        	continue
        if tstrand == "-":
        	qstart = qsize - qalg - qstart
        t = get_name(t)
        qend, tend = qstart + qalg, tstart + talg
        # get score, identity & overlap # LASTal is using +1/-1 for match/mismatch, while I need +1/0
        # calculate identity & overlap
        #identity = 1.0 * (score+(qalg-score)/2) / qalg
        #overlap  = 1.0 * qalg / qsize
        # store target & query info
        if t not in targets:
        	targets[t] = tsize 
        if q not in q2matches:
        	q2matches[q] = []
        	queries[q] = qsize
        # store
        q2matches[q].append((qstart, qend, t, tstart, tend, qstrand))
        
    # process queries by descending size
    tblocks = []
    print "# query\tbest target\taligned\tbest target size\t[%]\tquery aligned\tquery size\t[%]\tno of blocks"
    for i, q in enumerate(sorted(queries, key=lambda x: queries[x], reverse=True), 1): # 
	    # get best chromosome match
		bestt, talg, beststrand, total_alg = get_best_chr(q2matches[q])
		if bestt != 'CP002688': continue 
		blocks = get_blocks(q2matches[q], bestt, beststrand)
		print q, bestt, talg, targets[bestt], 100.*talg/targets[bestt], total_alg, queries[q], 100.*total_alg/queries[q], len(blocks)
		tblocks += blocks
		#'''
		for i, b in enumerate(blocks, 1): 
			print i, len(b), b[:2], b[-2:]
		return 
		#'''
    print "%s blocks for %s contigs"%(len(tblocks), len(queries))

def tab2errors():
	out = sys.stdout
	fasta, ref = sys.argv[1:3]
	process_hits(ref, fasta)    
            
if __name__=="__main__":
    t0 = datetime.now()
    tab2errors()
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)


