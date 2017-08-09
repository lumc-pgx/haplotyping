include: "helper.snake"
PARAMS = HaplotypingHelper(config, "Haplotype matching")

onsuccess: PARAMS.onsuccess()
onerror: PARAMS.onerror()

# main workflow
localrules:
    all

rule all:
    input:
        PARAMS.outputs

include: "rules/link_sources.snake"
include: "rules/matches.snake"
include: "rules/tabulate.snake"
include: "rules/pick.snake"
