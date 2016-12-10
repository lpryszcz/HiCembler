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

def intra_translocation_end(beststrand, tstart, blocks, lasti=-2):#+blocks[lasti][-1][2]
	if beststrand == "+" and tstart > blocks[lasti][-1][1] \
	  or beststrand == "-" and tstart < blocks[lasti][-1][1]: 
		return True
	return False

events = ["correct", "inversion", "intra-chr translocation", "inter-chr translocation"]
def get_blocks(matches, bestt, beststrand, maxtranfrac=0.1):
	"""Parse matches and report synteny blocks"""
	blocks, breaks = [], []
	# populate until bestt found
	for i, (qstart, qend, t, tstart, tend, strand) in enumerate(matches, 1):
		if t!=bestt: 
			if not blocks:
				breaks.append(events[3])
				blocks.append([])
			elif t!=blocks[-1][-1][0]:
				breaks.append(events[3])
				blocks.append([])
			blocks[-1].append((t, tstart, tend-tstart, strand))
		elif strand!=beststrand: 
			breaks.append(events[1])
			break
		else: 
			breaks.append(events[0])
			break
	blocks.append([(t, tstart, tend-tstart, strand)])

	lasti = -2
	large_intra_translocation = False
	for i in xrange(i, len(matches)):
		qstart, qend, t, tstart, tend, strand = matches[i]
		#print bestt, filter(lambda x: x[-1][0]==bestt, blocks)[-1][-1][-1], blocks
		#if i>5: break
		#if not filter(lambda x: x[-1][0]==bestt, blocks): print bestt, blocks
		# get end of intra-chr translocation
		if breaks[-1] == events[2]: 
			if intra_translocation_end(beststrand, tstart, blocks, lasti) or large_intra_translocation:
				# mark previous as translocation if smaller
				if sum(b[2] for b in blocks[lasti]) < sum(b[2] for b in blocks[-1]):
					#print "translocation changing order"
					breaks[lasti] = events[2]
					breaks[-1] = events[0]
				large_intra_translocation = False
				lasti = -2
			else:
				blocks[-1].append((t, tstart, tend-tstart, strand))
				continue
		# 3 inter-chromosomal translocation - this has to be first
		if t != matches[i-1][2]:
			blocks.append([])
			if t != bestt:
				breaks.append(events[3])
			elif strand != beststrand:
				breaks.append(events[1])				
			else:
				breaks.append(events[0])
		# 1 inversion - this has to be last; skip inter-chr translocations
		elif strand != filter(lambda x: x[-1][0]==bestt, blocks)[-1][-1][-1]:#matches[i-1][-1]:
			blocks.append([])
			if strand != beststrand: 
				breaks.append(events[1])
			else: 
				# reverse inverted block HERE!
				breaks.append(events[0])
		# 2 or intra-chromosomal translocation: tstart < previos for + or tstart > previous for -
		elif beststrand == "+" and tstart < filter(lambda x: x[-1][0]==bestt, blocks)[-1][-1][1] \
		  or beststrand == "-" and tstart > filter(lambda x: x[-1][0]==bestt, blocks)[-1][-1][1]:
		  	# split previous block, but ignore iter-chr translocation blocks
		  	lasti = [ii for ii, e in enumerate(breaks) if e!=events[3]][-1]
		  	if len(blocks[lasti])>1:
		  		lastblock = blocks[lasti]
		  		idx = [ii for ii, b in enumerate(lastblock, 1) \
		  				if beststrand=="+" and b[1]<tstart \
		  				or beststrand=="-" and b[1]>tstart]
		  		if idx and idx[-1] and idx[-1]<len(lastblock):
		  			blocks[lasti] = lastblock[:idx[-1]]
		  			blocks.append(lastblock[idx[-1]:])
		  			breaks.append(events[0])
			blocks.append([])
			breaks.append(events[2])
			# check for large intra-translocation ie. involving chr arm >10% of chr
			transsize = 0
			for m in filter(lambda x: x[2]==bestt, matches[i+1:]): 
				if intra_translocation_end(beststrand, m[3], blocks, lasti): 
					break
				transsize += m[4]-m[3]
			#if bestt=="CP002688": print transsize, sum(m[4]-m[3] for m in filter(lambda x: x[2]==bestt, matches))
			if transsize > maxtranfrac*sum(m[4]-m[3] for m in filter(lambda x: x[2]==bestt, matches)):
				large_intra_translocation = True
		# stop large intra-chr translocation with correct block
		elif breaks[-1] == events[2]:
			blocks.append([])
			breaks.append(events[0])
		'''	
		if blocks[-1]:
			# get gap size difference
			estgap = abs(qstart-matches[-1][1])
			trugap = abs(tstart-matches[-1][4])
		'''
		# add match
		blocks[-1].append((t, tstart, tend-tstart, strand))						
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
		if bestt != 'CP002684': continue 
		for i, (block, t) in enumerate(zip(blocks, breaks), 1): 
			print i, t, len(block), block[:2], block[-2:]
		print types, sum(types.values())
		print talg, sum(typesizes.values()), typesizes
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
        # skip human patches
        if t.startswith(("GL","KI")): continue
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


