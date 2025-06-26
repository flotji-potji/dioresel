rule picmin_genes:
    input:
        picmin = rules.picmin_to_bed.output,
        ann = "data/reference/vieillardii.gff"
    output:
        "raw_data/gene_func_ana/picmin_pcadapt_genes.bed"
    shell:
        r"""
        bedtools intersect -a {input.picmin} -b {input.ann} -wao \
            | awk '$11 ~ "gene"' \
            | cut -f1,2,3,4,5,6,7,8,17 \
            | awk -v OFS='\t' '{{split($9, a, "="); $9 = a[2]; print}}' \
            > {output}
        """

rule filter_genes:
    input:
        rules.picmin_genes.output
    output:
        "raw_data/gene_func_ana/filtered_genes.bed"
    shell:
        r"""
        awk '$8 < 0.05' {input} > {output}
        """

rule annotate_picmin:
    input:
        emapper = "data/emapper/emapper.annotations",
        picmin = rules.filter_genes.output
    output:
        "raw_data/gene_func_ana/picmin_pcadapt_genes.annotated.bed"
    params:
        col = 9
    shell:
        r"""
        awk -F'\t' -v OFS='\t' 'FNR==NR{{split($1,b,"."); 
                                        a[b[1]]; c[b[1]] = $8; d[b[1]] = $10}}; 
                                FNR!=NR{{gene = ${params.col};
                                        if(gene in a) 
                                        print $0, c[gene], d[gene]}}' \
            <(sort -k1 -V {input.emapper}) \
            <(sort -k{params.col} -V {input.picmin}) \
            > {output}
        """

rule shared_genes:
    input:
        rules.annotate_picmin.output
    output:
        "raw_data/gene_func_ana/shared_genes.bed"
    shell:
        r"""
        awk '$5 == 5' {input} > {output}
        """

rule fr_clade_genes:
    input:
        rules.annotate_picmin.output
    output:
        "raw_data/gene_func_ana/fr_genes.bed"
    shell:
        r"""
        awk '$5 <= 3' {input} > {output}
        """

rule sr_clade_genes:
    input:
        umb = "raw_data/mk_test/pair_cal_umb/mkt_table.tsv",
        vie = "raw_data/mk_test/pair_cal_vie/mkt_table.tsv"
    output:
        "raw_data/gene_func_ana/sr_genes.bed"
    shell:
        r"""
        cat {input.umb} {input.vie} \
            | awk -F'\t' '$7 > $6 && $8 < 0.05' \
            > {output}
        """

rule plot_gene_venn:
    input:
        rules.shared_genes.output,
        rules.fr_clade_genes.output,
        rules.sr_clade_genes.output
    output:
        "results/gene_func_ana/gene_venn.jpg"
    script:
        "../scripts/plot_gene_venn.R"

use rule annotate_picmin as annotate_sr_genes with:
    input:
        emapper = "data/emapper/emapper.annotations",
        picmin = rules.sr_clade_genes.output
    output:
        "raw_data/gene_func_ana/sr_genes.annotated.bed"
    params:
        col = 1

rule go_enrichment:
    input:
        shared = rules.shared_genes.output,
        fr_genes = rules.fr_clade_genes.output,
        sr_genes = rules.annotate_sr_genes.output
    output:
        "raw_data/gene_func_ana/go_enrichment.Robject",
        "raw_data/gene_func_ana/go_enrichment.tsv"
    script:
        "../scripts/topgo_enrichment.R"
