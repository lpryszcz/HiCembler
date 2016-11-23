### Table of Contents
- **[HiCembler](#HiCembler)**  
  - **[Installation](#Installation)**  
  - **[Running the pipeline](#running-the-pipeline)**  
    - **[Parameters](#parameters)**  
    - **[Test run](#test-run)**  
  - **[Support](#support)**
  - **[Citation](#citation)**
  
# HiCembler

This project was inspired by [DNA triangulation](https://github.com/NoamKaplan/dna-triangulation). 

## Installation
On most Linux distros, the installation should be as easy as:
```
sudo -H pip install -U matplotlib numpy scipy fastcluster pysam
git clone --recursive https://github.com/lpryszcz/redundans.git
cd HiCembler
(cd bin/snap && make clean && make)
```

### Dependencies
- `numpy`, `scipy`, `fastcluster`, `pysam` and `matplotlib` ie. `sudo -H pip install -U matplotlib numpy scipy fastcluster pysam`
- [SNAP aligner](https://github.com/amplab/snap)
- [sinkhorn_knopp](https://github.com/btaba/sinkhorn_knopp)


## Running the pipeline

### Parameters

### Test run

## Support 

## Citation
