rule intersect_all_pi_outliers:
    input:
        expand("raw_data/filtered_bed/pair_{p[0]}_{p[1]}/pi.outliers.bed", p=samples) 
    output:
        "raw_data/pi_within_intersect_tree/all_pi_outlier_intersect.bed" 
    shell:
        """
        bedtools multiinter -i {input} \
            | awk '$4 == 5' \
            | cut -f 1,2,3 \
            > {output}
        """

use rule intersect_all_pi_outliers as intersect_all_pair_pi_outliers with:
    input:
        expand("raw_data/filtered_bed/pair_{p[0]}_{p[1]}/pi_pair.outliers.bed", p=samples) 
    output:
        "raw_data/pi_between_intersect_tree/all_pi_outlier_intersect.bed" 

use rule pi_window_intersect as tree_pi_window_intersect with:
    input:
        pi = rules.intersect_all_pi_outliers.output,
        tajimad_within = rules.outliers_sel_tajimad.output,
        tajimad_between = rules.outliers_sel_pair_tajimad.output,
        pcadapt = rules.outliers_pcadapt.output,
        fst = rules.outliers_fst.output,
        mkt = rules.outliers_sel_mk_test.output
    output:
        "raw_data/pi_within_intersect_tree/pair_{sp1}_{sp2}/pi.tajimad.pcadapt.fst.mkt.bed"

use rule pi_window_intersect as tree_pair_pi_window_intersect with:
    input:
        pi = rules.intersect_all_pair_pi_outliers.output,
        tajimad_within = rules.outliers_sel_tajimad.output,
        tajimad_between = rules.outliers_sel_pair_tajimad.output,
        pcadapt = rules.outliers_pcadapt.output,
        fst = rules.outliers_fst.output,
        mkt = rules.outliers_sel_mk_test.output
    output:
        "raw_data/pi_between_intersect_tree/pair_{sp1}_{sp2}/pi.tajimad.pcadapt.fst.mkt.bed"

use rule pi_window_count as tree_pi_window_count with:
    input:
        rules.tree_pi_window_intersect.output
    output:
        "raw_data/pi_within_intersect_tree/pair_{sp1}_{sp2}/pi_summary.tsv"
    params:
        col = 4

use rule tree_pi_window_count as tree_pair_pi_window_count with:
    input:
        rules.tree_pair_pi_window_intersect.output
    output:
        "raw_data/pi_between_intersect_tree/pair_{sp1}_{sp2}/pi_summary.tsv"

use rule all_pi_counts as all_tree_pi_counts with:
    input:
        expand("raw_data/pi_within_intersect_tree/pair_{p[0]}_{p[1]}/pi_summary.tsv", p=samples) 
    output:
        "results/pi_within_intersect_tree/pi_summary.tsv"

use rule all_pi_counts as all_tree_pair_pi_counts with:
    input:
        expand("raw_data/pi_between_intersect_tree/pair_{p[0]}_{p[1]}/pi_summary.tsv", p=samples) 
    output:
        "results/pi_between_intersect_tree/pi_summary.tsv"

use rule plot_pi_intersect as plot_tree_pi_intersect with:
    input:
        rules.all_tree_pi_counts.output
    output:
        "results/pi_within_intersect_tree/pi_summary.jpg"
    params:
        plot_title = "shared PI-within windows"

use rule plot_pi_intersect as plot_tree_pair_pi_intersect with:
    input:
        rules.all_tree_pair_pi_counts.output
    output:
        "results/pi_between_intersect_tree/pi_summary.jpg"
    params:
        plot_title = "shared PI-between windows"