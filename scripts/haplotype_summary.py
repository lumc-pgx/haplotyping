import json

def tabulate_allele_matches(allele):
    """
    Convert allele match information into a tabular representation for printing
    """
    overview = "\n".join(
        [
            "- {}".format(allele["sequence_id"]),
            "- Variants: {}    Known: {}    Novel: {}    Significant: {}".format(
                len(allele["known_variants"]) + len(allele["novel_variants"]),
                len(allele["known_variants"]),
                len(allele["novel_variants"]),
                len(allele["significant_variants"])
            )
        ]
    )
    
    hits = allele["haplotypes"]
    columns = []
    columns.append(["{:25}".format("Variant")] + \
                   ["{:25}".format(v) for v in allele["known_variants"]] + \
                   ["{:25}".format(label) for label in ["Fraction", "Proportion", "Jaccard", "Significant"]])

    for hit in hits:
        column = ["{:6}".format(hit["haplotype"])]
        column += ["{:6}".format("$" if v in allele["significant_variants"] and v in hit["haplotype_variants"] else ("*" if v in hit["haplotype_variants"] else "")) for v in allele["known_variants"]]
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

