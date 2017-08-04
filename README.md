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

```plantuml
digraph snakemake_dag {
    graph[bgcolor=white, margin=0];
    node[shape=box, style=rounded, fontname=sans,                 fontsize=10, penwidth=2];
    edge[penwidth=2, color=grey];
	0[label = "all", color = "0.00 0.6 0.85", style="rounded"];
	1[label = "pick", color = "0.17 0.6 0.85", style="rounded"];
	2[label = "matches", color = "0.33 0.6 0.85", style="rounded"];
	3[label = "tabulate", color = "0.50 0.6 0.85", style="rounded"];
	3 -> 0
	1 -> 0
	2 -> 0
	2 -> 1
	2 -> 3
}
```
     
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
  - `pipe-runner`
- To specify that the pipeline should write output to a location other than the default:
  - `pipe-runner --directory path/to/output/directory`
        
