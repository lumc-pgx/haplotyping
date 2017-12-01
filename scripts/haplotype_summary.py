import json

def tabulate_allele_matches(allele):
    """
    Convert allele match information into a tabular representation for printing
    """
    g_notation = lambda variant: variant["g_notation"]
    is_significant = lambda variant: "significant" in variant["tags"]
    is_known = lambda variant: "known" in variant["tags"]
    is_novel = lambda variant: "novel" in variant["tags"]

    variants = allele["variants"]
    known_variants = [v for v in variants if is_known(v)]
    novel_variants = [v for v in variants if is_novel(v)]
    significant_variants = [v for v in variants if is_significant(v)]

    overview = "\n".join(
        [
            "- {}".format(allele["sequence_id"]),
            "- Variants: {}    Known: {}    Novel: {}    Significant: {}".format(
                len(variants),
                len(known_variants),
                len(novel_variants),
                len(significant_variants)
            )
        ]
    )
    
    hits = allele["haplotypes"]
    columns = []
    columns.append(["{:25}".format("Variant")] + \
                   ["{:25}".format(g_notation(v)) for v in known_variants] + \
                   ["{:25}".format(label) for label in ["Fraction", "Proportion", "Jaccard", "Significant"]])

    for hit in hits:
        column = ["{:6}".format(hit["haplotype"])]
        column += ["{:6}".format("$" if is_significant(v) and g_notation(v) in hit["haplotype_variants"] else \
                                     ("*" if g_notation(v) in hit["haplotype_variants"] else "")) for v in known_variants]
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


with open(snakemake.input.matches, 'r') as infile, open(snakemake.output[0], "w") as outfile:
    alleles = json.load(infile)
    results = ["\n".join(tabulate_allele_matches(a)) for a in alleles] 
    print("\n\n".join(results), file=outfile)

