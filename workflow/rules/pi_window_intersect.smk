rule filter_fst:
    input:
        rules.vcftools_to_bed.output
    output:
        "raw_data/filtered_bed/pair_{sp1}_{sp2}/fst.filtered.bed"
    shell:
        r"""
        awk '($5 > 3 && $4 > 0){{print}}' {input} > {output}
        """

rule outliers_fst:
    input:
        rules.filter_fst.output
    output:
        "raw_data/filtered_bed/pair_{sp1}_{sp2}/fst.outliers.bed"
    params:
        prob = 0.95,
        col = 4,
        sign = ">"
    shell:
        r"""
        quantile=$(cut -f {params.col} {input} \
            | Rscript -e 'print(quantile(scan("stdin", quiet=TRUE), {params.prob}))' \
            | tail -n +2)
        awk -v quant=$quantile '${params.col} {params.sign} quant {{print}}' {input} > {output}
        """

rule filter_pi:
    input:
        rules.pi_to_bed.output
    output:
        "raw_data/filtered_bed/pair_{sp1}_{sp2}/pi.filtered.bed"
    shell:
        r"""
        awk '($5 > 3 && $4 != "nan"){{print}}' {input} > {output}
        """

use rule filter_pi as filter_pair_pi with:
    input:
        rules.pair_pi_to_bed.output
    output:
        "raw_data/filtered_bed/pair_{sp1}_{sp2}/pi_pair.filtered.bed"

use rule outliers_fst as outliers_pi with:
    input:
        rules.filter_pi.output
    output:
        "raw_data/filtered_bed/pair_{sp1}_{sp2}/pi.outliers.bed"

use rule outliers_pi as outliers_pair_pi with:
    input:
        rules.filter_pair_pi.output
    output:
        "raw_data/filtered_bed/pair_{sp1}_{sp2}/pi_pair.outliers.bed"

use rule filter_pi as filter_tajimad with:
    input:
        rules.tajimad_to_bed.output
    output:
        "raw_data/filtered_bed/pair_{sp1}_{sp2}/tajimad.filtered.bed"

use rule filter_tajimad as filter_pair_tajimad with:
    input:
        rules.pair_tajimad_to_bed.output
    output:
        "raw_data/filtered_bed/pair_{sp1}_{sp2}/tajimad_pair.filtered.bed"

use rule outliers_fst as outliers_sel_tajimad with:
    input:
        rules.filter_tajimad.output
    output:
        "raw_data/filtered_bed/pair_{sp1}_{sp2}/tajimad.sel_outliers.bed"
    params:
        prob = 0.05,
        sign = "<"
        
use rule outliers_sel_tajimad as outliers_sel_pair_tajimad with:
    input:
        rules.filter_pair_tajimad.output
    output:
        "raw_data/filtered_bed/pair_{sp1}_{sp2}/tajimad_pair.sel_outliers.bed"

rule filter_pcadapt:
    input:
        rules.pcadapt_make_windows.output
    output:
        "raw_data/filtered_bed/pair_{sp1}_{sp2}/pcadapt.filtered.bed"
    shell:
        r"""
        awk '($6 > 3 && $5 != "NA" && $4 != "NA"){{print}}' {input} > {output}
        """

rule outliers_pcadapt:
    input:
        rules.filter_pcadapt.output
    output:
        "raw_data/filtered_bed/pair_{sp1}_{sp2}/pcadapt.outliers.bed"
    params:
        prob = 0.05
    shell:
        r"""
        awk '$4 < {params.prob} {{print}}' {input} > {output}
        """

rule outliers_sel_mk_test:
    input:
        rules.mk_test_to_window_bed.output
    output:
        "raw_data/filtered_bed/pair_{sp1}_{sp2}/mk_test.outliers.bed"
    params:
        prob = 0.05
    shell:
        r"""
        awk '($6 < {params.prob} && $5 > $4) {{print}}' {input} > {output}
        """

rule pi_window_intersect:
    input:
        pi = rules.outliers_pi.output,
        tajimad = rules.outliers_sel_tajimad.output,
        pcadapt = rules.outliers_pcadapt.output,
        fst = rules.outliers_fst.output,
        mkt = rules.outliers_sel_mk_test.output
    output:
        "raw_data/pi_within_intersect/pair_{sp1}_{sp2}/pi.tajimad.pcadapt.fst.mkt.bed"
    shell:
        r"""
        bedtools intersect -a {input.pi} \
                           -b {input.tajimad} {input.pcadapt} {input.fst} {input.mkt} \
                           -names tajimad pcadapt fst mkt \
                           -wo > {output}
        """

use rule pi_window_intersect as pair_pi_window_intersect with:
    input:
        pi = rules.outliers_pair_pi.output,
        tajimad = rules.outliers_sel_pair_tajimad.output,
        pcadapt = rules.outliers_pcadapt.output,
        fst = rules.outliers_fst.output,
        mkt = rules.outliers_sel_mk_test.output
    output:
        "raw_data/pi_between_intersect/pair_{sp1}_{sp2}/pi_pair.tajimad_pair.pcadapt.fst.mkt.bed"

rule pi_window_count:
    input:
        rules.pi_window_intersect.output
    output:
        "raw_data/pi_within_intersect/pair_{sp1}_{sp2}/pi_summary.tsv"
    shell:
        r"""
        cut -f 6 {input} \
            | sort \
            | bedtools groupby -g 1 -c 1 -o count \
            | awk -v OFS='\t' '$2 = $2 OFS "pair_{wildcards.sp1}_{wildcards.sp2}"' \
            > {output}
        """

use rule pi_window_count as pair_pi_window_count with:
    input:
        rules.pair_pi_window_intersect.output
    output:
        "raw_data/pi_between_intersect/pair_{sp1}_{sp2}/pi_summary.tsv"

rule all_pi_counts:
    input:
        expand("raw_data/pi_within_intersect/pair_{p[0]}_{p[1]}/pi_summary.tsv", p=samples) 
    output:
        "results/pi_within_intersect/pi_summary.tsv"
    shell:
        """
        cat {input} > {output}
        """

use rule all_pi_counts as all_pair_pi_counts with:
    input:
        expand("raw_data/pi_between_intersect/pair_{p[0]}_{p[1]}/pi_summary.tsv", p=samples) 
    output:
        "results/pi_between_intersect/pi_summary.tsv"

rule plot_pi_intersect:
    input:
        rules.all_pi_counts.output
    output:
        "results/pi_within_intersect/pi_summary.jpg"
    params:
        plot_title = "PI - within"
    script:
        "../scripts/plot_pi_intersect.R"

use rule plot_pi_intersect as plot_pair_pi_intersect with:
    input:
        rules.all_pair_pi_counts.output
    output:
        "results/pi_between_intersect/pi_summary.jpg"
    params:
        plot_title = "PI - between"
