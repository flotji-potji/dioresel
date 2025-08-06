rule pcadapt_intersection_upset:
    input:
       fg_genes = rules.go_enrichment_pcadapt_intersection.input.fg_genes,
       gen_table = rules.go_enrichment_pcadapt_intersection.output.table
    output:
        gene = "results/superexact/pcadapt_intersection/gene.pdf",
        func = "results/superexact/pcadapt_intersection/func.pdf",
    params:
        col = 1
    script:
        "../scripts/gene_functional_sharedness.R"

use rule pcadapt_intersection_upset as pcadapt_upset with:
    input:
       fg_genes = rules.go_enrichment_pcadapt.input.fg_genes,
       gen_table = rules.go_enrichment_pcadapt.output.table
    output:
        gene = "results/superexact/pcadapt/gene.pdf",
        func = "results/superexact/pcadapt/func.pdf",

use rule pcadapt_intersection_upset as pcadapt_mktest_upset with:
    input:
       fg_genes = rules.go_enrichment_pcadapt_mktest.input.fg_genes,
       gen_table = rules.go_enrichment_pcadapt_mktest.output.table
    output:
        gene = "results/superexact/pcadapt_mktest/gene.pdf",
        func = "results/superexact/pcadapt_mktest/func.pdf",