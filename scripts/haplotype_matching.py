import json
  
# make a set containing all reference variants for the gene
gene = snakemake.params.gene
reference_variants = {s["g_notation"] for s in gene["snps"]}
       
def characterise_allele(found_variants, significant_variants, gene_definition):
    """
    Compare the variants for an allele with all the known haplotypes.
    Characterise the commonality between the allele and the reference haplotypes.
    
    :param found_variants: A set of variants observed in the allele sequence
    :type found_variants: {str, str, ...} where str is a string containing an hgvs variant description
    
    :param significant_variants: A set of variants for the gene which are annotated as being 'significant'
    :type significant_variants: {str, str, ...} where str is a string containing an hgvs variant description

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
        haplotype_variants = {s["g_notation"] for s in gene_definition["snps"] if s["id"] in set(haplotype["snps"])}
        
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
                    "significant": list(shared_variants & significant_variants)
                }
            )
        except ZeroDivisionError:
            pass
    
    return summary


def match_allele(allele, gene, trim_boundary=False):
    """
    Compare the allele variants against the haplotype definitions
    Return a dictionary summarizing the comparison
    """
    # sort the variants
    variants = sorted(allele["variants"])
    
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
    significant_variants = {s["g_notation"] for s in gene["snps"] if s["tags"] is not None and \
        "significant" in s["tags"]} & found_variants

    # get all haplotype matches for this allele
    characterised_allele = characterise_allele(found_variants, significant_variants, gene)
    
    # filter to remove matches with zero overlap
    filtered_allele = [x for x in characterised_allele if x["fraction"] > 0]
    
    # sort by number of significant variants, fraction and jaccard
    # this gives a list of haplotypes ordered by decreasing 'score'
    sorted_allele = sorted(filtered_allele, key=lambda x: (len(x["significant"]), x["fraction"], x["jaccard"]), reverse=True)
    
    return {
        "sequence_id": allele["sequence_id"],
        "known_variants": sorted(list(known_variants)),
        "novel_variants": sorted(list(novel_variants)),
        "significant_variants": sorted(list(significant_variants)),
        "haplotypes":     sorted_allele
    }


with open(snakemake.input[0], 'r') as infile, open(snakemake.output[0], "w") as outfile:
    alleles = json.load(infile)
    results = [match_allele(allele, gene, snakemake.params.ignore_boundary) for allele in alleles]
    print(json.dumps(results, indent=4, separators=(',', ': ')), file=outfile)

