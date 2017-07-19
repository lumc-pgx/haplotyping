# load config file
configfile: srcdir("config.yaml")

# imports
import os
import glob
import yaml
import datetime

INPUT_FILES = glob.glob(config["VARIANT_DATA_PATH"] + "/*.json")
BARCODE_IDS = [".".join(os.path.basename(f).split(".")[:-1]) for f in INPUT_FILES]

with open(config["GENE"], "r") as infile:
    GENE = yaml.safe_load(infile)

# handlers for workflow exit status
onsuccess:
    print("Haplotype  workflow completed successfully")
    config_file = "config.{}.json".format("{:%Y-%m-%d_%H:%M:%S}".format(datetime.datetime.now()))
    with open(config_file, "w") as outfile:
        print(json.dumps(config), file=outfile)

onerror:
    print("Error encountered while executing workflow")
    shell("cat {log}")


# main workflow
localrules:
    all

rule all:
    input:
        expand("matches/{barcodes}.matches.json", barcodes=BARCODE_IDS),
        expand("tables/{barcodes}.matches.txt", barcodes=BARCODE_IDS),
        expand("haplotypes/{barcodes}.haplotype.txt", barcodes=BARCODE_IDS),
        "config.{}.yaml".format("{:%Y-%m-%d_%H:%M:%S}".format(datetime.datetime.now()))


rule matches:
    input:
        config["VARIANT_DATA_PATH"] + "/{barcode}.json"
    output:
        "matches/{barcode}.matches.json",
    params:
        gene = GENE,
        ignore_boundary = "OPTIONS" in config and "ignoreBoundary" in config["OPTIONS"]
    script:
        "scripts/haplotype_matching.py"


rule tabulate:
    input:
       "matches/{barcode}.matches.json",
    output:
       "tables/{barcode}.matches.txt"
    params:
        gene = GENE
    script:
        "scripts/haplotype_summary.py"


rule pick:
    input:
        "matches/{barcode}.matches.json",
    output:
        "haplotypes/{barcode}.haplotype.txt"
    params:
        default = GENE["haplotypes"][0]["type"]
    script:
        "scripts/pick_haplotype.py"


rule write_config:
    output:
        "config.{timestamp}.yaml"
    run:
        with open(output[0], "w") as outfile:
            yaml.dump(config, outfile, default_flow_style=False)
