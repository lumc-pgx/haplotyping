import json
from math import isclose
import locus_processing
  
# load the gene definition
locus = locus_processing.load_locus_yaml(snakemake.input.locus)
default = locus.haplotypes[0].type

def pick(allele):
    assignment = []
    if len(allele["haplotypes"]) > 0:
        top_jaccard = allele["haplotypes"][0]["jaccard"]
        top_significant = len(allele["haplotypes"][0]["significant"])

        assignment = [h["haplotype"] for h in allele["haplotypes"] if h["fraction"] == 1 and isclose(h["jaccard"], top_jaccard) and len(h["significant"]) == top_significant]
    
    if len(assignment) == 0:
        assignment = [default]

    return assignment


with open(snakemake.input.matches, 'r') as infile:
    alleles = json.load(infile)

for allele in alleles:
    allele["haplotype"] = pick(allele)

with open(snakemake.output[0], "w") as outfile:
    print(json.dumps(alleles, indent=4, separators=(',', ': ')), file=outfile)
