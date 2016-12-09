#!/usr/bin/env python
# Align contigs onto chromosomes, parse the alignments (.tab.gz) and report synteny blocks. 
# The difference in gap size will be calculated, but only within each synteny block. 

import os, sys, gzip
from datetime import datetime
#from FastaIndex import FastaIndex
from collections import Counter
import numpy as np
	
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

events = ["correct", "inversion", "intra-chr translocation", "inter-chr translocation"]
def get_blocks(matches, bestt, beststrand):
	"""Parse matches and report synteny blocks"""
	qstart, qend, t, tstart, tend, strand = matches[0]
	blocks = [[(t, tstart, tend-tstart, strand),]]#matches[0][2:6],]]
	breaks = [events[0]]
	event = False
	reverse = False
	for i in xrange(1, len(matches)):
		qstart, qend, t, tstart, tend, strand = matches[i]#; print matches[i]
		# 1 inversion 
		if strand != matches[i-1][-1]:
			blocks.append([])
			if strand!=beststrand: 
				breaks.append(events[1])
			else:
				breaks.append(events[0])
		# 2 or intra-chromosomal translocation: tstart < previos for + or tstart > previous for -
		elif strand == "+" and tstart < matches[i-1][3] or strand == "-" and tstart > matches[i-1][3]:
			blocks.append([])
			breaks.append(events[2])
			# need to get end of intra-translocation! but how?
		# 3 inter-chromosomal translocation
		elif t != matches[i-1][2]:
			blocks.append([])
			if t!=bestt:
				breaks.append(events[3])
			else:
				breaks.append(events[0])
		elif blocks[-1]:
			# get gap size difference
			estgap = abs(qstart-matches[-1][1])
			trugap = abs(tstart-matches[-1][4])
		# add match
		blocks[-1].append((t, tstart, tend-tstart, strand))#matches[i][2:6])
	return blocks, breaks

def process_hits(ref, fasta):
	"""Report rearrengements from assembly"""
	q2matches, queries, targets = get_matches(ref, fasta)
    # process queries by descending size
	tblocks = []
	header = "# query\tbest target\taligned\tbest target size\t[%]\tquery aligned\tquery size\t[%]\t# blocks"
	for e in events: 
		header += "\t# %s\t%s [bp]\t[%s]" % (e, e, '%')
	print header
	totdata = np.zeros(header.count('\t')-1)
	for i, q in enumerate(sorted(queries, key=lambda x: queries[x], reverse=True), 1): # 
	    # get best chromosome match
		bestt, talg, beststrand, total_alg = get_best_chr(q2matches[q])
		# get synteny blocks and event types
		blocks, breaks = get_blocks(q2matches[q], bestt, beststrand)
		# get no. of each event and its size
		tblocks += blocks
		types = {e: 0 for e in events} #Counter(breaks)
		typesizes = {e: 0 for e in events}
		for i, (block, t) in enumerate(zip(blocks, breaks), 1): 
			types[t] += 1
			typesizes[t] += sum(b[2] for b in block)
		'''
		if bestt != 'CP002688': continue 
		for i, (block, t) in enumerate(zip(blocks, breaks), 1): 
			print i, len(block), block[:2], block[-2:]
		print types, sum(types.values())
		print talg, sum(typesizes.values()), typesizes
		return 
		#'''
		# report
		data = [q, bestt, talg, targets[bestt], 100.*talg/targets[bestt]]
		data += [total_alg, queries[q], 100.*total_alg/queries[q], len(blocks)]
		for e in events:
			data += [types[e], typesizes[e], 100.*typesizes[e]/sum(typesizes.values())]
		print "\t".join(map(str, data))
		totdata += data[2:]
	# recaulculate %
	totdata[2] = 100.*totdata[0]/totdata[1]
	totdata[5] = 100.*totdata[3]/totdata[4]
	indices = [8,11,14,17]#; print sum(totdata[indices])
	for i in indices: 
		totdata[i+1] = 100.*totdata[i]/sum(totdata[indices])
	print "SUM\t-\t"+"\t".join(map(str, totdata))

def get_name(name): return name.split()[0].split('|')[-1].split('.')[0]

def get_matches(ref, fasta, identityTh=0.9, overlapTh=.66, minscore=100, threads=4):
    """Generate & process LASTal hits"""
    # execute last and get best query-to-reference matches only
    tabfn = generate_tab(ref, fasta, threads)
    q2matches, queries, targets = {}, {}, {}
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
    return q2matches, queries, targets
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


