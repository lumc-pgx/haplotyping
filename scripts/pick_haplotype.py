import json
from math import isclose
 
def pick(allele):
    assignment = []
    if len(allele["haplotypes"]) > 0:
        top_jaccard = allele["haplotypes"][0]["jaccard"]
        top_significant = len(allele["haplotypes"][0]["significant"])

        assignment = [h["haplotype"] for h in allele["haplotypes"] if h["fraction"] == 1 and isclose(h["jaccard"], top_jaccard) and len(h["significant"]) == top_significant]
    
    if len(assignment) == 0:
        assignment = [snakemake.params.default]

    return assignment


with open(snakemake.input[0], 'r') as infile, open(snakemake.output[0], "w") as outfile:
    alleles = json.load(infile)
    results = ["\t".join([a["sequence_id"]] + pick(a)) for a in alleles]
    print("\n".join(results), file=outfile)

