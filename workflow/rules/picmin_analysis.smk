rule intersect_pcadapt_fst:
    input:
        pcadapt = rules.pcadapt_make_windows.output,
        fst = rules.vcftools_to_bed.output
    output:
        pcadapt = "raw_data/picmin/pair_{sp1}_{sp2}/pcadapt.intersect.bed",
        fst = "raw_data/picmin/pair_{sp1}_{sp2}/fst.intersect.bed"
    shell:
        r"""
        bedtools intersect -a {input.pcadapt} -b {input.fst} -wb > {output.pcadapt}
        cut -f 7,8,9,10,11 {output.pcadapt} > {output.fst}
        """

rule scale_pcadapt:
    input:
        fst = rules.filter_fst.output,
        pcadapt = rules.filter_pcadapt.output
    output:
        "raw_data/picmin/pair_{sp1}_{sp2}/pcadapt.scaled.bed"
    shell:
        r"""
        fst_mean=$(awk '{{total=$4}} END{{print total/NR}}' {input.fst}) 
        awk -v mean=$fst_mean -v OFS='\t' '$5 = $5/mean' {input.pcadapt} > {output}
        """

use rule filter_pcadapt as filter_scaled_pcadapt with:
    input:
        rules.scale_pcadapt.output
    output:
        "raw_data/picmin/pair_{sp1}_{sp2}/pcadapt.scaled.filtered.bed"

rule picmin_reduce_bed_and_create_df_pcadapt:
    input:
        rules.sort_pcadapt.output
    output:
        "raw_data/picmin/pair_{sp1}_{sp2}/pcadapt_picmin.Robject"
    params:
        pval_small = True
    script:
        "../scripts/picmin_create_df.R"

use rule picmin_reduce_bed_and_create_df_pcadapt as picmin_reduce_bed_and_create_df_fst with:
    input:
        rules.filter_fst.output
    output:
        "raw_data/picmin/pair_{sp1}_{sp2}/fst_picmin.Robject"

use rule picmin_reduce_bed_and_create_df_pcadapt as picmin_reduce_bed_and_create_df_tajimad with:
    input:
        rules.filter_tajimad.output
    output:
        "raw_data/picmin/pair_{sp1}_{sp2}/tajimad_picmin.Robject"
    params:
        pval_small = False

rule picmin_analysis_pcadapt:
    input:
        expand("raw_data/picmin/pair_{p[0]}_{p[1]}/pcadapt_picmin.Robject", p=samples)
    output:
        "raw_data/picmin/picmin_pcadapt_results.Robject"
    params:
        num_reps = 100000
    script:
        "../scripts/picmin_analysis.R"

rule picmin_to_bed:
    input:
        rules.picmin_analysis_pcadapt.output
    output:
        "raw_data/picmin/picmin_pcadapt_results.bed"
    script:
        "../scripts/picmin_to_bed.R"

use rule picmin_analysis_pcadapt as picmin_analysis_fst with:
    input:
        expand("raw_data/picmin/pair_{p[0]}_{p[1]}/fst_picmin.Robject", p=samples)
    output:
        "raw_data/picmin/picmin_fst_results.Robject"

use rule picmin_analysis_pcadapt as picmin_analysis_tajimad with:
    input:
        expand("raw_data/picmin/pair_{p[0]}_{p[1]}/tajimad_picmin.Robject", p=samples)
    output:
        "raw_data/picmin/picmin_tajimad_results.Robject"

rule plot_picmin_pcadapt_results:
    input:
        rules.picmin_analysis_pcadapt.output
    output:
        "results/picmin/pcadapt/picmin_results.pval.jpeg",
        "results/picmin/pcadapt/picmin_results.padj.jpeg",
        "results/picmin/pcadapt/picmin_results.lin_venn.jpeg",
        "results/picmin/pcadapt/picmin_results.est_venn.jpeg",
        "results/picmin/pcadapt/picmin_results.lin_est_heat.jpeg"
    script:
        "../scripts/picmin_plot_results.R"

use rule plot_picmin_pcadapt_results as plot_picmin_fst_results with:
    input:
        rules.picmin_analysis_fst.output
    output:
        "results/picmin/fst/picmin_results.pval.jpeg",
        "results/picmin/fst/picmin_results.padj.jpeg",
        "results/picmin/fst/picmin_results.lin_venn.jpeg",
        "results/picmin/fst/picmin_results.est_venn.jpeg",
        "results/picmin/fst/picmin_results.lin_est_heat.jpeg"

use rule plot_picmin_pcadapt_results as plot_picmin_tajimad_results with:
    input:
        rules.picmin_analysis_tajimad.output
    output:
        "results/picmin/tajimad/picmin_results.pval.jpeg",
        "results/picmin/tajimad/picmin_results.padj.jpeg",
        "results/picmin/tajimad/picmin_results.lin_venn.jpeg",
        "results/picmin/tajimad/picmin_results.est_venn.jpeg",
        "results/picmin/tajimad/picmin_results.lin_est_heat.jpeg"
