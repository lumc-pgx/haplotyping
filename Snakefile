# load config file
configfile: srcdir("config.yaml")

# imports
import os
import glob
import yaml

INPUT_FILES = glob.glob(config["VARIANT_DATA_PATH"] + "/*.json")
BARCODE_IDS = [".".join(os.path.basename(f).split(".")[:-1]) for f in INPUT_FILES]

with open(config["GENE"], "r") as infile:
    GENE = yaml.safe_load(infile)

# handlers for workflow exit status
onsuccess:
    print("Haplotype  workflow completed successfully")
onerror:
    print("Error encountered while executing workflow")
    shell("cat {log}")


# main workflow
localrules:
    all

rule all:
    input:
        expand("haplotypes/{barcodes}.haplotypes.json", barcodes=BARCODE_IDS)


rule haplotypes:
    input:
        alleles = config["VARIANT_DATA_PATH"] + "/{barcode}.json"
    output:
        json = "haplotypes/{barcode}.haplotypes.json",
        table = "haplotypes/{barcode}.haplotypes.txt"
    params:
        gene = GENE,
        ignore_boundary = "OPTIONS" in config and "ignoreBoundary" in config["OPTIONS"]
    script:
        "scripts/compare_haplotypes.py"



