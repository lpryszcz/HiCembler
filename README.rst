.. contents:: Table of Contents

HiCembler
=========


This project was inspired by `DNA triangulation <https://github.com/NoamKaplan/dna-triangulation>`_. 

============
Installation
============
On most Linux distros, the installation should be as easy as:


.. code-block:: bash

    sudo -H pip install -U matplotlib numpy scipy fastcluster pysam ete3
    git clone --recursive https://github.com/lpryszcz/redundans.git
    cd HiCembler
    (cd bin/snap && make clean && make)
    (cd bin/idba && ./build.sh && ./configure && make)


Dependencies
~~~~~~~~~~~~
- Python dependencies: matplotlib, numpy, scipy, sklearn, fastcluster, pysam, ete3, `sinkhorn_knopp <https://github.com/btaba/sinkhorn_knopp>`_
.. code-block:: bash
                  
   sudo -H pip install -U matplotlib numpy scipy sklearn fastcluster pysam ete3 sinkhorn_knopp


- `SNAP aligner <https://github.com/amplab/snap>`_
- `IDBA <https://github.com/loneknightpy/idba>`_ optionally

====================
Running the pipeline
====================




Parameters
~~~~~~~~~~


  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -v, --verbose         verbose
  -i BAM, --bam BAM
                        BAM file(s)
  -f FASTA, --fasta FASTA
                        Genome FastA file
  -o OUTDIR, --outdir OUTDIR
                        output name
  -w WINDOWSIZE, --windowSize WINDOWSIZE
                        window size in kb used for scaffolding [[5, 10, 2]]
  -m MINSIZE, --minSize MINSIZE
                        minimum contig length [2000]
  -q MAPQ, --mapq MAPQ  mapping quality [10]
  -u UPTO, --upto UPTO  process up to this number of reads from each library [all]
  -t THREADS, --threads THREADS
                        number of processes to use [4]
  -d DPI, --dpi DPI     output images dpi [300]
  --minWindows MINWINDOWS
                        minimum number of windows per contig, has to be >1 [2]


Test run
~~~~~~~~



=======
Support
=======



========
Citation
========
