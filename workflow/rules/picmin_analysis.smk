rule picmin_reduce_bed_and_create_df_pcadapt:
    input:
        rules.filter_pcadapt.output
    output:
        "raw_data/picmin/pair_{sp1}_{sp2}/pcadapt_picmin.Robject"
    params:
        pval_small = True
    script:
        "../scripts/picmin_create_df.R"

rule picmin_analysis_pcadapt:
    input:
        expand("raw_data/picmin/pair_{p[0]}_{p[1]}/pcadapt_picmin.Robject", p=samples)
    output:
        "raw_data/picmin/picmin_pcadapt_results.Robject"
    params:
        num_reps = 100000
    script:
        "../scripts/picmin_analysis.R"

rule plot_picmin_pcadapt_results:
    input:
        rules.picmin_analysis_pcadapt.output
    output:
        "raw_data/picmin/picmin_pcadapt_results.pval.jpeg",
        "raw_data/picmin/picmin_pcadapt_results.padj.jpeg"
    script:
        "../scripts/picmin_plot_manhattan.R"
