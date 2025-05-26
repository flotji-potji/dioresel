rule pi_window_intersect:
    input:
        pi = rules.outliers_pi.output,
        tajimad_within = rules.outliers_sel_tajimad.output,
        tajimad_between = rules.outliers_sel_pair_tajimad.output,
        pcadapt = rules.outliers_pcadapt.output,
        fst = rules.outliers_fst.output,
        mkt = rules.outliers_sel_mk_test.output
    output:
        "raw_data/pi_within_intersect/pair_{sp1}_{sp2}/pi.tajimad.pcadapt.fst.mkt.bed"
    shell:
        r"""
        bedtools intersect -a {input.pi} \
                           -b {input.tajimad_within} \
                              {input.tajimad_between} {input.pcadapt} \
                              {input.fst} {input.mkt} \
                           -names tajimad_within tajimad_between pcadapt fst mkt \
                           -wo > {output}
        """

use rule pi_window_intersect as pi_bottom_window_intersect with:
    input:
        pi = rules.outliers_pi_bottom.output,
        tajimad_within = rules.outliers_sel_tajimad.output,
        tajimad_between = rules.outliers_sel_pair_tajimad.output,
        pcadapt = rules.outliers_pcadapt.output,
        fst = rules.outliers_fst.output,
        mkt = rules.outliers_sel_mk_test.output
    output:
        "raw_data/pi_within_intersect/pair_{sp1}_{sp2}/pi_bottom.tajimad.pcadapt.fst.mkt.bed"

use rule pi_window_intersect as pair_pi_window_intersect with:
    input:
        pi = rules.outliers_pair_pi.output,
        tajimad_within = rules.outliers_sel_tajimad.output,
        tajimad_between = rules.outliers_sel_pair_tajimad.output,
        pcadapt = rules.outliers_pcadapt.output,
        fst = rules.outliers_fst.output,
        mkt = rules.outliers_sel_mk_test.output
    output:
        "raw_data/pi_between_intersect/pair_{sp1}_{sp2}/pi_pair.tajimad_pair.pcadapt.fst.mkt.bed"

use rule pi_window_intersect as pair_pi_bottom_window_intersect with:
    input:
        pi = rules.outliers_pair_pi_bottom.output,
        tajimad_within = rules.outliers_sel_tajimad.output,
        tajimad_between = rules.outliers_sel_pair_tajimad.output,
        pcadapt = rules.outliers_pcadapt.output,
        fst = rules.outliers_fst.output,
        mkt = rules.outliers_sel_mk_test.output
    output:
        "raw_data/pi_between_intersect/pair_{sp1}_{sp2}/pi_pair_bottom.tajimad_pair.pcadapt.fst.mkt.bed"

rule pi_window_count:
    input:
        rules.pi_window_intersect.output
    output:
        "raw_data/pi_within_intersect/pair_{sp1}_{sp2}/pi_summary.tsv"
    params:
        col = 6
    shell:
        r"""
        cut -f {params.col} {input} \
            | sort \
            | bedtools groupby -g 1 -c 1 -o count \
            | awk -v OFS='\t' '$2 = $2 OFS "pair_{wildcards.sp1}_{wildcards.sp2}"' \
            > {output}
        """

use rule pi_window_count as pi_bottom_window_count with:
    input:
        rules.pi_bottom_window_intersect.output
    output:
        "raw_data/pi_within_intersect/pair_{sp1}_{sp2}/pi_bottom_summary.tsv"

use rule pi_window_count as pair_pi_window_count with:
    input:
        rules.pair_pi_window_intersect.output
    output:
        "raw_data/pi_between_intersect/pair_{sp1}_{sp2}/pi_summary.tsv"

use rule pi_window_count as pair_pi_bottom_window_count with:
    input:
        rules.pair_pi_bottom_window_intersect.output
    output:
        "raw_data/pi_between_intersect/pair_{sp1}_{sp2}/pi_bottom_summary.tsv"

rule all_pi_counts:
    input:
        expand("raw_data/pi_within_intersect/pair_{p[0]}_{p[1]}/pi_summary.tsv", p=samples) 
    output:
        "results/pi_within_intersect/pi_summary.tsv"
    shell:
        """
        cat {input} > {output}
        """

use rule all_pi_counts as all_pi_bottom_counts with:
    input:
        expand("raw_data/pi_within_intersect/pair_{p[0]}_{p[1]}/pi_bottom_summary.tsv", p=samples) 
    output:
        "results/pi_within_intersect/pi_bottom_summary.tsv"


use rule all_pi_counts as all_pair_pi_counts with:
    input:
        expand("raw_data/pi_between_intersect/pair_{p[0]}_{p[1]}/pi_summary.tsv", p=samples) 
    output:
        "results/pi_between_intersect/pi_summary.tsv"

use rule all_pi_counts as all_pair_pi_bottom_counts with:
    input:
        expand("raw_data/pi_between_intersect/pair_{p[0]}_{p[1]}/pi_bottom_summary.tsv", p=samples) 
    output:
        "results/pi_between_intersect/pi_bottom_summary.tsv"

rule plot_pi_intersect:
    input:
        rules.all_pi_counts.output
    output:
        "results/pi_within_intersect/pi_summary.jpg"
    params:
        plot_title = "PI - within"
    script:
        "../scripts/plot_pi_intersect.R"

use rule plot_pi_intersect as plot_pi_bottom_intersect with:
    input:
        rules.all_pi_bottom_counts.output
    output:
        "results/pi_within_intersect/pi_bottom_summary.jpg"
    params:
        plot_title = "PI - bottom windows - within"

use rule plot_pi_intersect as plot_pair_pi_intersect with:
    input:
        rules.all_pair_pi_counts.output
    output:
        "results/pi_between_intersect/pi_summary.jpg"
    params:
        plot_title = "PI - between"

use rule plot_pi_intersect as plot_pair_pi_bottom_intersect with:
    input:
        rules.all_pair_pi_bottom_counts.output
    output:
        "results/pi_between_intersect/pi_bottom_summary.jpg"
    params:
        plot_title = "PI - bottom windows - between"
