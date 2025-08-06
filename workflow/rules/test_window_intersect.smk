rule goi_intersect_all:
    # goi = genes of interest
    input:
        goi = "data/paun_et_al/genes_of_interest.bed",
        pcadapt = rules.outliers_pcadapt.output,
        mktest = rules.outliers_sel_mk_test.output,
        fst = rules.outliers_fst.output,
        pi = rules.outliers_pi.output,
        pi_bottom = rules.outliers_pi_bottom.output,
        pair_pi = rules.outliers_pair_pi.output,
        pair_pi_bottom = rules.outliers_pair_pi_bottom.output,
        tajimad = rules.outliers_sel_tajimad.output,
        pair_tajimad = rules.outliers_sel_pair_tajimad.output,
    output:
        "raw_data/test_intersect/pair_{sp1}_{sp2}/goi_intersect.bed"
    params:
        names = "pi pi_bottom pair_pi pair_pi_bottom tajimad pair_tajimad fst pcadapt mktest"
    shell:
        r"""
        bedtools intersect -a {input.goi} \
                           -b {input.pi} {input.pi_bottom} \
                              {input.pair_pi} {input.pair_pi_bottom} \
                              {input.tajimad} {input.pair_tajimad} \
                              {input.fst} \
                              {input.pcadapt} {input.mktest} \
                           -names {params.names} -wo \
            | cut -f1,2,3,4,5,6 \
            > {output}
        """        

use rule pi_window_count as goi_window_count with:
    input:
        rules.goi_intersect_all.output
    output:
        "raw_data/test_intersect/pair_{sp1}_{sp2}/goi_count.tsv"
    params:
        col = 6

use rule all_pi_counts as all_goi_counts with:
    input:
        expand("raw_data/test_intersect/pair_{p[0]}_{p[1]}/goi_count.tsv", p=samples) 
    output:
        "results/test_intersect/goi_count.merged.tsv"
    
rule plot_goi_counts:
    input:
        rules.all_goi_counts.output
    output:
        "results/test_intersect/goi_intersection_heatmap.svg"
    script:
        "../scripts/plot_goi_heatmap.R"
    
rule intersect_all_tests:
    input:
        annotation = "results/gene_annotation/func_gene_annotation.windows.bed",
        pcadapt = rules.outliers_pcadapt.output,
        mktest = rules.outliers_sel_mk_test.output,
        fst = rules.outliers_fst.output,
        pi = rules.outliers_pi.output,
        pi_bottom = rules.outliers_pi_bottom.output,
        pair_pi = rules.outliers_pair_pi.output,
        pair_pi_bottom = rules.outliers_pair_pi_bottom.output,
        tajimad = rules.outliers_sel_tajimad.output,
        pair_tajimad = rules.outliers_sel_pair_tajimad.output,
    output:
        "raw_data/test_intersect/pair_{sp1}_{sp2}/all_tests_intersect.bed"
    shell:
        r"""
        bedtools multiinter -i {input} \
            | bedtools expand -c 5 \
            | sort -k5,5n \
            | bedtools groupby -g 5 -c 6,7,8,9,10,11,12,13,14,15 \
                               -o sum,sum,sum,sum,sum,sum,sum,sum,sum,sum \
            > {output}
        """

rule plot_all_tests_intersection:
    input:
        expand("raw_data/test_intersect/pair_{p[0]}_{p[1]}/all_tests_intersect.bed", p=samples) 
    output:
        "results/test_intersect/all_tests_intersection.svg"
    script:
        "../scripts/plot_all_tests_heatmap.R"

rule pcadapt_intersect_all:
    input:
        pcadapt = rules.outliers_pcadapt.output,
        pi = rules.outliers_pi.output,
        pi_bottom = rules.outliers_pi_bottom.output,
        pair_pi = rules.outliers_pair_pi.output,
        pair_pi_bottom = rules.outliers_pair_pi_bottom.output,
        tajimad = rules.outliers_sel_tajimad.output,
        pair_tajimad = rules.outliers_sel_pair_tajimad.output,
        fst = rules.outliers_fst.output
    output:
        "raw_data/test_intersect/pair_{sp1}_{sp2}/pcadapt.intersect.bed"
    params:
        names = "pi pi_bottom pair_pi pair_pi_bottom tajimad pair_tajimad fst"
    shell:
        r"""
        bedtools intersect -a {input[0]} \
                           -b {input.pi} {input.pi_bottom} \
                              {input.pair_pi} {input.pair_pi_bottom} \
                              {input.tajimad} {input.pair_tajimad} \
                              {input.fst} \
                           -names {params.names} \
                           -wo > {output}
        """

use rule pcadapt_intersect_all as mktest_intersect_all with:
    input:
        mktest = rules.outliers_sel_mk_test.output,
        pi = rules.outliers_pi.output,
        pi_bottom = rules.outliers_pi_bottom.output,
        pair_pi = rules.outliers_pair_pi.output,
        pair_pi_bottom = rules.outliers_pair_pi_bottom.output,
        tajimad = rules.outliers_sel_tajimad.output,
        pair_tajimad = rules.outliers_sel_pair_tajimad.output,
        fst = rules.outliers_fst.output
    output:
        "raw_data/test_intersect/pair_{sp1}_{sp2}/mktest.intersect.bed"

use rule pi_window_count as pcadapt_window_count with:
    input:
        rules.pcadapt_intersect_all.output
    output:
        "raw_data/test_intersect/pair_{sp1}_{sp2}/pcadapt_count.tsv"
    params:
        col = 7

use rule pi_window_count as mktest_window_count with:
    input:
        rules.mktest_intersect_all.output
    output:
        "raw_data/test_intersect/pair_{sp1}_{sp2}/mktest_count.tsv"
    params:
        col = 8

rule pcadapt_window_number:
    input:
        pcadapt = rules.outliers_pcadapt.output,
        summary = rules.pcadapt_window_count.output
    output:
        "raw_data/test_intersect/pair_{sp1}_{sp2}/pcadapt_summary.tsv"
    params:
        test = "pcadapt"
    shell:
        r"""
        awk -v OFS='\t' 'FNR!=NR && FNR == 1 {{a=NR-1}} FNR!=NR {{print $0, "{params.test}", a}}' \
            {input} > {output}
        """

use rule pcadapt_window_number as mktest_window_number with:
    input:
        pcadapt = rules.outliers_sel_mk_test.output,
        summary = rules.mktest_window_count.output
    output:
        "raw_data/test_intersect/pair_{sp1}_{sp2}/mktest_summary.tsv"
    params:
        test = "mktest"

rule merge_test_ouputs:
    input:
        rules.sort_fst.output,
        rules.sort_pi.output,
        rules.sort_pair_pi.output,
        rules.sort_tajimad.output,
        rules.sort_pair_tajimad.output,
    output:
        "raw_data/filtered_bed/pair_{sp1}_{sp2}/tests.merged.tsv"
    shell:
        r"""
        i=0
        for file in {input}; do
            awk -v OFS='|' -v file=$file \
                'BEGIN{{split(file, a, "/"); split(a[4], b, "."); print b[1], b[1]"_sp"}}
                {{print $4, "pair_{wildcards.sp1}_{wildcards.sp2}"}}' \
                $file > $TMPDIR/tmp${{i}}
            i=$((i+1))
        done
        paste $TMPDIR/tmp* > {output}
        rm $TMPDIR/tmp*
        """

use rule all_pi_counts as all_merged_test_outputs with:
    input:
        expand("raw_data/filtered_bed/pair_{p[0]}_{p[1]}/tests.merged.tsv", p=samples) 
    output:
        "results/filtered_bed/tests.merged.tsv"

use rule all_pi_counts as all_pcadapt_summaries with:
    input:
        expand("raw_data/test_intersect/pair_{p[0]}_{p[1]}/pcadapt_summary.tsv", p=samples) 
    output:
        "results/test_intersect/pcadapt_summary.tsv"

use rule all_pi_counts as all_mktest_summaries with:
    input:
        expand("raw_data/test_intersect/pair_{p[0]}_{p[1]}/mktest_summary.tsv", p=samples) 
    output:
        "results/test_intersect/mktest_summary.tsv"

use rule all_pi_counts as merge_all_test_summaries with:
    input:
        rules.all_pcadapt_summaries.output,
        rules.all_mktest_summaries.output
    output:
        "results/test_intersect/test_summary.tsv"

rule plot_test_intersect:
    input:
        test_outputs = rules.all_merged_test_outputs.output,
        test_summary = rules.merge_all_test_summaries.output
    output:
        "results/test_intersect/pcadapt/fst.svg",
        "results/test_intersect/pcadapt/pi_between.svg",
        "results/test_intersect/pcadapt/pi_within.svg",
        "results/test_intersect/pcadapt/tajimad.svg"
    script:
        "../scripts/plot_test_intersect.R"
        


