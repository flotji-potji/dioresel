
rule vcftools_step_window:
    # rule produces sliding window results as output
    input:
        vcf="data/variants/pair_{sp1}_{sp2}/vieref_{sp1}_{sp2}.vcf.gz",
        sp1_sample="data/variants/pair_{sp1}_{sp2}/{sp1}.samples",
        sp2_sample="data/variants/pair_{sp1}_{sp2}/{sp2}.samples"
    output:
        "raw_data/fst_step_window/pair_{sp1}_{sp2}/vieref_{sp1}_{sp2}.windowed.weir.fst"
    params:
        output_format="raw_data/fst_step_window/pair_{sp1}_{sp2}/vieref_{sp1}_{sp2}",
        window_size=10000
    shell:
        r"""
        vcftools --gzvcf {input.vcf} \
		--weir-fst-pop {input.sp1_sample} \
		--weir-fst-pop {input.sp2_sample} \
		--fst-window-size {params.window_size} \
		--out {params.output_format}
        """

rule vcftools_to_bed:
    input:
        rules.vcftools_step_window.output
    output:
        "raw_data/fst_step_window/pair_{sp1}_{sp2}/vieref_{sp1}_{sp2}.bed"
    params:
        header = "WEIGHTED_FST",
        window_size = 10000,
        add_one = 0
    shell:
        r"""
        awk 'BEGIN{{col = 0; nv_col = 0; OFS = "\t"}} NR==1{{
                                    for (i=1;i<=NF;i++) {{
                                        if ($i == "{params.header}") {{
                                            col = i
                                        }}
                                        if ($i == "N_VARIANTS") {{
                                            nv_col = i
                                        }}
                                    }}
                                }} NR>1{{print $1, $2+{params.add_one}, $2+{params.window_size}+{params.add_one}, $col, $nv_col}}' \
            {input} > {output}
        """

rule r_plot_fst:
    input:
        rules.vcftools_step_window.output
    output:
        "results/fst_plots/pair_{sp1}_{sp2}/pair_{sp1}_{sp2}.fst.jpg"
    params:
        plot_title="Fst - {sp1} vs. {sp2}",
        window_info="window size: 10000bp"
    script:
        "../scripts/plot_fst.R"