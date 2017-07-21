# Haplotyping workflow

**Comparison of allelic variants and haplotype definitions**

The haplotyping workflow module performs the following operations:  
- Identification of allele variants shared with reference haplotypes.
- Prioritisation of haplotype matches.
- Identification of the highest scoring haplotype.

The pipeline outputs three files per barcode:
- A {barcode}.matches.json file which contains all of the match information per barcode
- A {barcode}.matches.txt file which contains a human readable summary of the per-allele matches
- A {barcode}.haplotype.txt file which contains the haplotype assignment per allele

## Requirements
- [Conda/Miniconda](https://conda.io/miniconda.html)  

## Installation
- Clone the repository
  - `git clone https://git.lumc.nl/PharmacogenomicsPipe/haplotyping.git`

- Change to the haplotyping directory
  - `cd haplotyping`

- Create a conda environment for running the pipeline
  - `conda env create -n haplotyping -f environment.yaml`

- In order to use the pipeline on the cluster, update your .profile to use the drmaa library:
  - `echo "export DRMAA_LIBRARY_PATH=libdrmaa.so.1.0" >> ~/.profile`
  - `source ~/.profile`

## Configuration
Pipeline configuration settings can be altered by editing [config.yaml](config.yaml).  

## Execution
- Activate the conda environment
  - `source activate haplotyping`
- For parallel execution on the cluster
  - `./run_cluster.sh`
- To specify that the pipeline should write output to a location other than the default:
  - `./run_cluster.sh -d path/to/output/directory`

