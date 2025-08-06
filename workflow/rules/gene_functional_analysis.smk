rule annotate_pcadapt_intersection:
    input:
        intersections = expand("raw_data/test_intersect/pair_{p[0]}_{p[1]}/pcadapt.intersect.bed", p=samples),
        annotation = rules.func_annotation_to_window_bed.output
    output:
        "raw_data/gene_func_ana/pcadapt_intersection.annotated.bed"
    params:
        spp_names = expand("pair_{p[0]}_{p[1]}", p=samples),
        test_col = ",17"
    shell:
        r"""
        bedtools intersect -a {input.annotation} \
                           -b {input.intersections} \
                           -names {params.spp_names} \
                           -wo \
            | cut -d$'\t' -f 1,2,3,4,5,6,7,8,9,10{params.test_col} \
            > {output}
        """

use rule annotate_pcadapt_intersection as annotate_pcadapt with:
    input:
        intersections = expand("raw_data/filtered_bed/pair_{p[0]}_{p[1]}/pcadapt.outliers.bed", p=samples),
        annotation = rules.func_annotation_to_window_bed.output
    output:
        "raw_data/gene_func_ana/pcadapt.annotated.bed"
    params:
        test_col = ""

use rule annotate_pcadapt_intersection as annotate_mktest with:
    input:
        intersections = expand("raw_data/filtered_bed/pair_{p[0]}_{p[1]}/mk_test.outliers.bed", p=samples),
        annotation = rules.func_annotation_to_window_bed.output
    output:
        "raw_data/gene_func_ana/mktest.annotated.bed"
    params:
        test_col = ""

rule annotate_picmin:
    input:
        intersections = rules.picmin_to_bed.output,
        annotation = rules.func_annotation_to_window_bed.output
    output:
        "raw_data/gene_func_ana/picmin.annotated.bed"
    shell:
        r"""
        bedtools intersect -a {input.annotation} \
                           -b <(awk '$7 < 0.05' {input.intersections}) \
                           -wo \
            | cut -d$'\t' -f 1,2,3,4,5,6,7,8,9,13 \
            > {output}
        """

rule merge_pcadapt_mktest_annotation:
    input:
        pcadapt = rules.annotate_pcadapt_intersection.output,
        mktest = rules.annotate_mktest.output,
    output:
        "raw_data/gene_func_ana/pcadapt_mktest.annotated.bed"
    shell:
        r"""
        grep fst {input.pcadapt} \
            | grep 'cal_spn\|cal_eru\|cal_ruf' \
            | cut -f1,2,3,4,5,6,7,8,9,10 \
            > {output}
        grep 'cal_umb\|cal_vie' {input.mktest} \
            >> {output}
        """

rule go_enrichment_pcadapt_intersection:
    input:
        fg_genes = rules.annotate_pcadapt_intersection.output,
        bg_genes = rules.merge_emapper_interpro.output
    output:
        fst = "results/gene_func_ana/pcadapt_intersection/fst_enrichment.jpg",
        pair_pi = "results/gene_func_ana/pcadapt_intersection/pair_pi_enrichment.jpg",
        pi = "results/gene_func_ana/pcadapt_intersection/pi_enrichment.jpg",
        pi_bottom = "results/gene_func_ana/pcadapt_intersection/pi_bottom_enrichment.jpg",
        tajimad = "results/gene_func_ana/pcadapt_intersection/tajimad_enrichment.jpg",
        table = "raw_data/gene_func_ana/pcadapt_intersection/enrichment_table.Robject"
    params:
        plot_title = "pcadapt intersected with"
    script:
        "../scripts/topgo_enrichment.R"

use rule go_enrichment_pcadapt_intersection as go_enrichment_pcadapt with:
    input:
        fg_genes = rules.annotate_pcadapt.output,
        bg_genes = rules.merge_emapper_interpro.output
    output:
        unfiltered = "results/gene_func_ana/pcadapt/enrichment.jpg",
        table = "raw_data/gene_func_ana/pcadapt/enrichment_table.Robject"
    params:
        plot_title = "pcadapt"

use rule go_enrichment_pcadapt_intersection as go_enrichment_mktest with:
    input:
        fg_genes = rules.annotate_mktest.output,
        bg_genes = rules.merge_emapper_interpro.output
    output:
        unfiltered = "results/gene_func_ana/mktest/enrichment.jpg",
        table = "raw_data/gene_func_ana/mktest/enrichment_table.Robject"
    params:
        plot_title = "MKT"

use rule go_enrichment_pcadapt_intersection as go_enrichment_pcadapt_mktest with:
    input:
        fg_genes = rules.merge_pcadapt_mktest_annotation.output,
        bg_genes = rules.merge_emapper_interpro.output
    output:
        unfiltered = "results/gene_func_ana/pcadapt_mktest/enrichment.jpg",
        table = "raw_data/gene_func_ana/pcadapt_mktest/enrichment_table.Robject"
    params:
        plot_title = "pcadapt/fst - MKT"

use rule go_enrichment_pcadapt_intersection as go_enrichment_picmin with:
    input:
        fg_genes = rules.annotate_picmin.output,
        bg_genes = rules.merge_emapper_interpro.output
    output:
        unfiltered = "results/gene_func_ana/picmin/enrichment.jpg",
        table = "raw_data/gene_func_ana/picmin/enrichment_table.Robject"
    params:
        plot_title = "picmin"
