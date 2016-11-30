#!/usr/bin/env python
desc="""Report contact matrix and contact plot. 
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Bratislava, 19/10/2016
"""

# Force matplotlib to not use any Xwindows backend.
import matplotlib; matplotlib.use('Agg')

import gzip, os, resource, sys
import commands, os, subprocess, sys
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
from matplotlib.ticker import FuncFormatter, MaxNLocator
from FastaIndex import FastaIndex
from multiprocessing import Pool
from collections import Counter

def logger(message, log=sys.stdout):
    """Log messages"""
    memory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    log.write("[%s] %s    [memory: %6i Mb]\n"%(datetime.ctime(datetime.now()), message, memory))

def fasta2windows(fasta, windowSize, verbose, skipShorter=1, minSize=2000,
                  contigs=[], check=1, filterwindows=1):
    """Generate windows over chromosomes"""
    # init fasta index
    faidx = FastaIndex(fasta)
    # 
    if verbose:
        logger("Parsing FastA file...")
    # filter windows so they are smaller than largest chr and withing reasonalbe range toward genome size
    if filterwindows:    
        maxchrlen = max(faidx.id2stats[c][0] for c in faidx) # 100-20000 windows
        windowSize = filter(lambda x: 1000*x<maxchrlen/10 and 1000*x<0.01*faidx.genomeSize and 1000*x>0.00005*faidx.genomeSize, windowSize)
        if verbose:
            logger(" selected %s windows [kb]: %s"%(len(windowSize), str(windowSize)))
    windowSize = [w*1000 for w in windowSize]
    # generate windows
    windows, chr2window = [[] for w in windowSize], [{} for w in windowSize]
    genomeSize = 0
    base2chr = {}
    skipped = []
    for i, c in enumerate(faidx, 1): 
        if contigs and c not in contigs:
            continue
        if i%1e5 == 1:
            sys.stderr.write(' %s   \r'%i)
        # get windows
        size = faidx.id2stats[c][0]
        # skip short contigs    
        if skipShorter and size < min(windowSize) or size<minSize:
            skipped.append(size)
            continue
        for i in range(len(windowSize)):
            # skip contig if shorter than given window
            if skipShorter and size < windowSize[i]:
                #print size, "skipped %s"%windowSize[i]
                continue
            # get starting window
            chr2window[i][c] = len(windows[i])
            for start in range(0, size, windowSize[i]):
                windows[i].append((c, start, start+windowSize[i]))
            # update last entry end
            if not skipShorter:
                windows[i][-1] = (c, start, size)
        # get chromosome tick in the middle    
        base2chr[genomeSize+size/2] = c
        # update genomeSize
        genomeSize += size
    # print [len(w) for w in windows]
    if verbose:
        logger(' %s bases in %s contigs divided in %s-%s windows. '%(faidx.genomeSize, len(faidx), len(windows[0]), len(windows[-1])))
        if skipped:
            logger('  %s bases in %s contigs skipped.'%(sum(skipped), len(skipped)))
    return windowSize, windows, chr2window, base2chr, faidx.genomeSize

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

def get_window(chrom, start, length, windowSize, chr2window):
    """Return window alignment belongs to. """
    end = start + length
    # get middle position
    pos = int(start+round((end-start)/2))
    # get window
    window = pos / windowSize + chr2window[chrom]
    return window

def _bam2array(args):
    """Parse pysam alignments and return matching windows"""
    _bam, regions, mapq, windowSize, chr2window, upto = args 
    data = [Counter() for ii in range(len(windowSize))]
    proc = _get_samtools_proc(_bam, mapq, regions)
    # process reads from given library
    for i, line in enumerate(proc.stdout, 1):
        if upto and i > upto:
            break
        # unload data
        ref1, start1, mapq, cigar, ref2, start2, insertSize, seq = line.split('\t')[2:10]
        start1, start2, seqlen = int(start1), int(start2), len(seq)
        # update ref if alignment in the same chrom
        if ref2 == "=":
            ref2 = ref1
        # skip if contig not present in array
        if ref1 not in chr2window[0] or ref2 not in chr2window[0]:
            continue
        # get windows
        for ii in range(len(windowSize)):
            w1 = get_window(ref1, start1, seqlen, windowSize[ii], chr2window[ii])
            w2 = get_window(ref2, start2, seqlen, windowSize[ii], chr2window[ii])
            # matrix is symmetrix, so make sure to add only to one part
            if w2 < w1:
                w1, w2 = w2, w1
            # update contact array
            data[ii][(w1, w2)] += 1
                
    # stop subprocess
    proc.terminate()
    return data

def bam2array_multi(windows, windowSize, chr2window, bam, mapq, upto=0,
                    regions=[], verbose=1, minSize=2000, threads=4):
    """Return contact matrix based on BAM"""
    arrays = [np.zeros((len(w), len(w)), dtype="float32") for w in windows]
    i = 1
    if not regions:
        regions = chr2window[-1].keys()
    if threads > len(regions):
        threads = len(regions)
    # init empty array
    p = Pool(threads)
    args = [(_bam, regions[ii::threads], mapq, windowSize, chr2window, upto) for ii in range(threads) for _bam in bam]
    for wdata in p.imap_unordered(_bam2array, args): 
        for ii in range(len(wdata)):
            for (w1, w2), c in wdata[ii].iteritems():
                arrays[ii][w1][w2] += c
                i += c
        if upto and i>upto:
            break
        sys.stderr.write(" %s \r"%i)
    if verbose:
        logger(" %s alignments parsed"%i)
    return arrays

def plot(outfn, a, genomeSize, base2chr, _windowSize, dpi=300, ext="svg"):
    """Save contact plot"""
    
    def format_fn(tick_val, tick_pos):
        """Mark axis ticks with chromosome names"""
        if int(tick_val) in base2chr:
            return base2chr[int(tick_val)]
        else:
            sys.stderr.write("[WARNING] %s not in ticks!\n"%tick_val)
            return ''
            
    # invert base2chr
    base2chr = {genomeSize-b: c for b, c in base2chr.iteritems()}
    # start figure
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title("Contact intensity plot [%sk]"%(_windowSize/1000,))
    # label Y axis with chromosome names
    if len(base2chr)<50:
        ax.yaxis.set_major_formatter(FuncFormatter(format_fn))
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        plt.yticks(base2chr.keys())
        ax.set_ylabel("Chromosomes")
    else:
        ax.set_ylabel("Genome position")
    # label axes
    ax.set_xlabel("Genome position")        
    plt.imshow(a+1, cmap=cm.hot, norm=LogNorm(), extent=(0, genomeSize, 0, genomeSize))# 
    plt.colorbar()
    # save
    fig.savefig("%s.%s"%(outfn,ext), dpi=dpi, papertype="a4")

def bam2plot(fasta, bam, outbase, windowSize, mapq=10, threads=1,
                upto=float('inf'), dpi=300, dtype='float32', verbose=1):
    """Convert SAM to SSPACE TAB file."""
    # get windows
    windowSize, windows, chr2window, base2chr, genomeSize = fasta2windows(fasta, windowSize, verbose)

    # get contact matrices
    logger("Parsing alignements...")
    arrays = bam2array_multi(windows, windowSize, chr2window, bam, mapq, upto=upto, threads=threads)
    
    # save
    logger("Saving & plotting...")
    for a, _windowSize, _windows in zip(arrays, windowSize, windows):
        _outfn = os.path.join(outbase, "%sk"%(_windowSize/1000,))
        if verbose:
            logger(" %s"%_outfn)
        # save windows
        with gzip.open(_outfn+".windows.tab.gz", "w") as out:
            out.write("\n".join("\t".join(map(str, w)) for w in _windows)+"\n")

        #save normalised
        with open(_outfn+".npz", "w") as out:
            np.savez_compressed(out, a)

        # make symmetric
        a += a.T - np.diag(a.diagonal())
            
        if len(_windows)<5e4:
            plot(_outfn, a, genomeSize, base2chr, _windowSize, dpi)
        elif verbose:
            sys.stderr.write("[WARNING] Very large matrix (%s x %s). Skipped plotting!\n"%(len(_windows), len(_windows)))
    logger("Done!")

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--version', action='version', version='1.0b')   
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose")    
    parser.add_argument("-i", "--bam", nargs="+", help="BAM files")
    parser.add_argument("-f", "--fasta", type=file, help="Genome FastA file")
    parser.add_argument("-o", "--output", required=1, help="output name")
    parser.add_argument("-w", "--windowSize", nargs="+", default=[10000, 1000, 500, 200, 100, 50, 20, 10, 5, 2, 1], type=int,
                        help="window size in kb [%(default)s]")
    parser.add_argument("-m", "--mapq", default=10, type=int,
                        help="mapping quality [%(default)s]")
    parser.add_argument("-u", "--upto", default=0, type=float,
                        help="process up to this number of reads from each library [all]")
    parser.add_argument("-t", "--threads", default=4, type=int,
                        help="number of processes to use [%(default)s]")
    parser.add_argument("-d", "--dpi", default=300, type=int,
                        help="output images dpi [%(default)s]")
    parser.add_argument("--dtype", default='float32', 
                        help="numpy array data type (try uint16 if MemoryError for 100+k windows) [%(default)s]")

    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
        
    # check if output not present
    fileexists = 0
    for w in o.windowSize:
        outfn = "%s.%sk.npz"%(o.output, w)
        if os.path.isfile(outfn):
            sys.stderr.write("Output file exists: %s !\n"%outfn)
            fileexists = 1
    if fileexists:
        sys.exit(1)

    # create outdir
    if not os.path.isdir(o.output):
        os.makedirs(o.output)
        
    # process
    bam2plot(o.fasta, o.bam, o.output, o.windowSize, o.mapq, o.threads,
             o.upto, o.dpi, o.dtype, o.verbose)

if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    #except IOError as e:
    #    sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror))
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
    