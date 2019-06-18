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
        graph [bb="0,0,140.6,252",
                bgcolor=white,
                margin=0
        ];
        node [fontname=sans,
                fontsize=10,
                label="\N",
                penwidth=2,
                shape=box,
                style=rounded
        ];
        edge [color=grey,
                penwidth=2
        ];
        0        [color="0.00 0.6 0.85",
                height=0.5,
                label=all,
                pos="40.6,18",
                width=0.75];
        1        [color="0.13 0.6 0.85",
                height=0.5,
                label=tabulate,
                pos="40.6,90",
                width=0.76389];
        1 -> 0   [pos="e,40.6,36.104 40.6,71.697 40.6,63.983 40.6,54.712 40.6,46.112"];
        2        [color="0.27 0.6 0.85",
                height=0.5,
                label=link_source,
                pos="40.6,234",
                width=0.97222];
        3        [color="0.40 0.6 0.85",
                height=0.5,
                label=matches,
                pos="40.6,162",
                width=0.8125];
        2 -> 3   [pos="e,40.6,180.1 40.6,215.7 40.6,207.98 40.6,198.71 40.6,190.11"];
        3 -> 0   [pos="e,24.564,36.154 24.564,143.85 16.587,134.12 7.7754,121.27 3.6004,108 -1.2001,92.737 -1.2001,87.263 3.6004,72 6.699,62.148 12.352,\
52.53 18.321,44.251"];
        3 -> 1   [pos="e,40.6,108.1 40.6,143.7 40.6,135.98 40.6,126.71 40.6,118.11"];
        4        [color="0.53 0.6 0.85",
                height=0.5,
                label=pick,
                pos="113.6,90",
                width=0.75];
        3 -> 4   [pos="e,95.763,108.1 58.645,143.7 67.662,135.05 78.718,124.45 88.544,115.03"];
        4 -> 0   [pos="e,58.438,36.104 95.555,71.697 86.539,63.05 75.482,52.449 65.657,43.027"];
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
        
