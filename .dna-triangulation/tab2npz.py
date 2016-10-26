#!/usr/bin/env python
# Convert tab / tab.gz to .npz

import gzip, sys
import numpy as np
from triangulation import _load_raw

file_txt = sys.argv[1]

if file_txt.endswith('.gz'):
    fh = gzip.open(file_txt, 'r')
    d, bin_chr, bin_position = _load_raw(fh)
else:
    fh = open(file_txt, 'r')
    d, bin_chr, bin_position = _load_raw(fh)

np.savez_compressed(file_txt+".npz", d)
with gzip.open(file_txt+".windows.tab.gz", "w") as out:
    for chrom, (start, end) in zip(bin_chr, bin_position):
        out.write("%s\t%s\t%s\n"%(chrom, start, end))