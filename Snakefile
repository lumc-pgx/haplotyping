include: "helper.snake"
PARAMS = Haplotyping(config, "Haplotype matching")

onsuccess: PARAMS.onsuccess()
onerror: PARAMS.onerror()

# main workflow
localrules:
    all


rule all:
    input:
        PARAMS.outputs


rule matches:
    input:
        config["VARIANT_DATA_PATH"] + "/{barcode}.json"
    output:
        "matches/{barcode}.matches.json",
    params:
        gene = PARAMS.GENE,
        ignore_boundary = PARAMS.IGNORE_BOUNDARY 
    script:
        "scripts/haplotype_matching.py"


rule tabulate:
    input:
       "matches/{barcode}.matches.json",
    output:
       "tables/{barcode}.matches.txt"
    params:
        gene = PARAMS.GENE
    script:
        "scripts/haplotype_summary.py"


rule pick:
    input:
        "matches/{barcode}.matches.json",
    output:
        "haplotypes/{barcode}.haplotype.txt"
    params:
        default = PARAMS.GENE["haplotypes"][0]["type"]
    script:
        "scripts/pick_haplotype.py"

