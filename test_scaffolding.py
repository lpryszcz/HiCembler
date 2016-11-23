
import os, sys
root = os.path.dirname(os.path.abspath(sys.argv[0]))
src = ["/home/lpryszcz/src/HiCembler/bin", "/home/lpryszcz/src/HiCembler/bin/snap", "/home/lpryszcz/src/HiCembler/bin/sinkhorn_knopp"]
paths = [os.path.join(root, p) for p in src]
sys.path = paths + sys.path
os.environ["PATH"] = "%s:%s"%(':'.join(paths), os.environ["PATH"])
from bam2scaffolds import *
from array2scaffolds import _average_reduce, _average_reduce_2d

def join_scaffolds(i, scaffold1, scaffold2, d, contig2indices, minWindows=3):
    """Join two adjacent scaffolds"""
    indices1 = get_indices(scaffold1, contig2indices)
    indices2 = get_indices(scaffold2, contig2indices)
    # skip contigs with less windows than minWindows
    if len(indices1) < len(indices2):
        scaffold1, indices1, scaffold2, indices2 = scaffold2, indices2, scaffold1, indices1
    if len(indices2) < minWindows:
        return scaffold1
    # get subset of array for scaffold1 and scaffold2
    _d = d[:, indices1+indices2][indices1+indices2, :]
    # get orientation: 0: s-s; 1: s-e; 2: e-s; 3: e-e
    ## contact values for ends of two contigs are compared
    ## and the max value is taken as true contact
    ### this should be accurate, but you may consider some ML function
    n1, n2 = len(indices1), len(indices2)
    # compare diagonals of contact matrix
    i = n2/2
    #dflip = np.fliplr(_d[:n1,n2:])
    #d1, d2, d3, d4 = np.diag(_d, k=n1), np.diag(dflip), np.diag(dflip, k=n2-n1), np.diag(_d, k=n2)
    orientation = np.argmax(map(sum, (d1[:i], d2[:i], d3[-i:], d4[-i:])))
    # s - s
    if   orientation == 0:
        scaffold = get_reversed(scaffold2) + scaffold1 #get_reversed(scaffold1) + scaffold2
    # s - e
    elif orientation == 1:
        scaffold = scaffold2 + scaffold1 # get_reversed(scaffold1) + get_reversed(scaffold2)
    # e - s
    elif orientation == 2:
        scaffold = scaffold1 + scaffold2 
    # e - e
    else:
        scaffold = scaffold1 + get_reversed(scaffold2)
    '''header = "%s x %s\torientation: %s\tshape: %s"%(n1, n2, orientation, _d.shape[0])
    footer = " ".join(n for n, rev in scaffold)
    fname = "tmp/%7s.tsv"%i; fname = fname.replace(' ','0')
    np.savetxt(fname, _d, fmt='%5.f', delimiter='\t', newline='\n', header=header, footer=footer, comments='# ')'''
    return scaffold

bam = ["/home/lpryszcz/cluster/hic/arath/platanus/ARATH.d05.l100_contig.fa.bam"]
fasta ='/home/lpryszcz/cluster/hic/arath/platanus/ARATH.d05.l100_contig.fa'
outdir = '/home/lpryszcz/cluster/hic/arath/platanus/bam2scaffolds.new_diag'
windowSize = [100]
windowSize2 = [10]
minWindows = 2
mapq = 0
dpi = 100
upto = 0
verbose = 1
i = 1

# get clusters
clusters, _windowSize = bam2clusters(bam, fasta, outdir, windowSize, mapq, dpi, upto, verbose); print [len(c) for c in clusters]

# scaffolding
contigs = clusters[0][i-1]

windowSize2, windows, chr2window, base2chr, genomeSize = fasta2windows(fasta, windowSize2, verbose=0, skipShorter=1, contigs=set(contigs), filterwindows=0)

arrays = [np.zeros((len(w), len(w)), dtype="float32") for w in windows]
arrays = bam2array(arrays, windowSize2, chr2window, bam,  mapq, regions=contigs, verbose=1)

d = arrays[0]
# generate missing handles
bin_chr = np.array([c for c, s, e in windows[0]])
bin_position = np.array([(s, e) for c, s, e in windows[0]])

# make symmetric & normalise
#d = normalize(d)
d += d.T; d -= np.diag(d.diagonal()/2); d = normalize_rows(d)
#d, bin_chr, bin_position = normalize_diagional(d, bin_chr, bin_position)
logger(" cluster_%s with %s windows in %s contigs"%(i, d.shape[0], len(contigs)))

# get tree on reduced matrix
t = array2tree(transform(_average_reduce_2d(d, bin_chr)), np.unique(bin_chr))
#t = array2tree(transform(d), bin_chr)

# get scaffold
#scaffold = tree2scaffold(t, d, bin_chr, bin_position, minWindows)

contig2indices = get_contig2indices(bin_chr)

# populate internal nodes with growing scaffolds
for ii, n in enumerate(t.traverse('postorder'), 1):
    # add scaffold for each leave
    # each scaffold consists of ordered list of contigs and their orientations (0/False: Fwd or 1/True:Rev)
    if n.is_leaf():
        n.scaffold = [(n.name, 0)]
        continue
    # unload children
    n1, n2 = n.get_children()
    # and combine scaffolds
    n.scaffold = join_scaffolds(i, n1.scaffold, n2.scaffold, d, contig2indices, minWindows)
    #if ii>10: break
    
--

