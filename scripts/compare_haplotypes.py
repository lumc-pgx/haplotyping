import json, yaml
  
# make a set containing all reference variants for the gene
gene = snakemake.params.gene
reference_variants = {gene["snps"][snp]["g_notation"] for snp in gene["snps"]}
       
def characterise_allele(found_variants, gene_definition):
    """
    Compare the variants for an allele with all the known haplotypes.
    Characterise the commonality between the allele and the reference haplotypes.
    
    :param found_variants: A set of variants observed in the allele sequence
    :type found_variants: {str, str, ...} where str is a string containing an hgvs variant description
    
    :param gene_definition: The parsed gene file containing snp and haplotype definitions
    :type gene_definition: dict
    
    :return: a list of dicts where each dict characterises the commonality between the allelic
             variants and the variants for a single haplotype
    """
    # list to hold the comparisons 
    summary = []
    
    # compare the allele variants with those for each defined haplotype
    for haplotype in (h for h in gene_definition["haplotypes"] if len(h["snps"]) > 0):
        # the variants for this haplotype 
        haplotype_variants = {gene_definition["snps"][s]["g_notation"] for s in haplotype["snps"]}
        
        # the variants shared between the allele and the haplotype
        shared_variants = found_variants & haplotype_variants
        
        # summarize
        try:
            summary.append(
                {
                    # the current haplotype
                    "haplotype": haplotype["type"], 
                    
                    # the variants for this haplotype
                    "haplotype_variants": sorted(list(haplotype_variants)), 
                    
                    # the variants shared between allele and haplotype
                    "shared_variants": sorted(list(shared_variants)), 
                    
                    # the fraction of haplotype variants observed
                    "fraction": len(shared_variants)/ len(haplotype_variants),
                    
                    # jaccard coefficient of haplotype and allele variants
                    "jaccard": len(shared_variants) / (len(haplotype_variants) + len(found_variants) - len(shared_variants)),
                   
                    # shared 'significant' variants
                    "significant": [
                        gene_definition["snps"][s]["g_notation"] for s in haplotype["snps"] if gene_definition["snps"][s]["tags"] is not None and \
                        "significant" in gene_definition["snps"][s]["tags"] and \
                        gene_definition["snps"][s]["g_notation"] in shared_variants
                    ]
                }
            )
        except ZeroDivisionError:
            pass
    
    return summary


def match_allele(variants, gene, trim_boundary=False):
    """
    Compare the allele variants against the haplotype definitions
    Return a dictionary summarizing the comparison
    """
    # sort the variants
    variants = sorted(variants)
    
    # trim the first and last variants if requested
    if trim_boundary:
        if len(variants) > 0 and "ins" in variants[0]:
            variants = variants[1:]
    
        if len(variants) > 0  and "ins" in variants[-1]:
            variants = variants[:-1]
    
    # set of variants found in the allele
    found_variants = set(variants)
    
    # known variants - intersection of found variants with the set of all reference variants
    known_variants = found_variants & reference_variants
    
    # novel variants found in allele but not found in the reference set
    novel_variants = found_variants - reference_variants
    
    # significant variants found in allele
    significant_variants = {
        gene["snps"][s]["g_notation"] for s in gene["snps"] if gene["snps"][s]["tags"] is not None and \
        "significant" in gene["snps"][s]["tags"] and \
        gene["snps"][s]["g_notation"] in found_variants
    }

    # get all haplotype matches for this allele
    characterised_allele = characterise_allele(found_variants, gene)
    
    # filter to remove matches with zero overlap
    filtered_allele = [x for x in characterised_allele if x["fraction"] > 0]
    
    # sort by number of significant variants, fraction and jaccard
    # this gives a list of haplotypes ordered by decreasing 'score'
    sorted_allele = sorted(filtered_allele, key=lambda x: (len(x["significant"]), x["fraction"], x["jaccard"]), reverse=True)
    
    return {
        "known_variants": sorted(list(known_variants)),
        "novel_variants": sorted(list(novel_variants)),
        "significant_variants": sorted(list(significant_variants)),
        "haplotypes":     sorted_allele
    }


def tabulate_allele_matches(allele, match):
    """
    Convert allele match information into a tabular representation for printing
    """
    overview = "\n".join(
        [
            "- {}    %id: {}".format(allele["sequence_id"], allele["identity"]),
            "- Variants: {}    Known: {}    Novel: {}    Significant: {}".format(
                len(match["known_variants"]) + len(match["novel_variants"]),
                len(match["known_variants"]),
                len(match["novel_variants"]),
                len(match["significant_variants"])
            )
        ]
    )
    
    hits = match["haplotypes"]
    columns = []
    columns.append(["{:25}".format("Variant")] + \
                   ["{:25}".format(v) for v in match["known_variants"]] + \
                   ["{:25}".format(label) for label in ["Fraction", "Proportion", "Jaccard", "Significant"]])

    for hit in hits:
        column = ["{:6}".format(hit["haplotype"])]
        column += ["{:6}".format("$" if v in match["significant_variants"] and v in hit["haplotype_variants"] else ("*" if v in hit["haplotype_variants"] else "")) for v in match["known_variants"]]
        column.append("{:6}".format("{}/{}".format(len(hit["shared_variants"]), len(hit["haplotype_variants"]))))
        column.append("{:6}".format(str(hit["fraction"])[:4]))
        column.append("{:6}".format(str(hit["jaccard"])[:4]))
        column.append("{:<6}".format(len(hit["significant"])))
        columns.append(column)
        
    results = []
    for row in range(len(columns[0])):
        results.append("".join([c[row] for c in columns]))
    
    table = "\n".join(results)
    
    return (overview, table)


with open(snakemake.input.alleles, 'r') as infile:
    alleles = json.load(infile)

results = []
for allele in alleles:
    match = match_allele(allele["variants"], gene, snakemake.params.ignore_boundary)
    match["sequence_id"] = allele["sequence_id"]
    
    table = tabulate_allele_matches(allele, match)

    assignment = []
    if len(match["haplotypes"]) > 0:
        top_jaccard = match["haplotypes"][0]["jaccard"]
        top_significant = match["haplotypes"][0]["significant"]

        assignment = [h["haplotype"] for h in match["haplotypes"] if h["fraction"] == 1 and h["jaccard"] == top_jaccard and h["significant"] == top_significant]
    
    if len(assignment) == 0:
        assignment = [gene["haplotypes"][0]["type"]]

    results.append(
        (
            match,
            "\n".join(table),
            "\t".join([allele["sequence_id"]] + assignment)
        )
    )

with open(snakemake.output.json, "w") as outfile:
    print(json.dumps([r[0] for r in results], indent=4, separators=(',', ': ')), file=outfile)

with open(snakemake.output.table, "w") as outfile:
    print("\n\n".join([r[1] for r in results]), file=outfile)

with open(snakemake.output.call, "w") as outfile:
    print("\n".join([r[2] for r in results]), file=outfile)
