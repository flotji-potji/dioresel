
rule vcftools_step_window:
    # rule produces sliding window results as output
    input:
        vcf="data/variants/pair_{sp1}_{sp2}/vieref_{sp1}_{sp2}.vcf.gz",
        sp1_sample="data/variants/pair_{sp1}_{sp2}/{sp1}.samples",
        sp2_sample="data/variants/pair_{sp1}_{sp2}/{sp2}.samples"
    output:
        "raw_data/fst_step_window/pair_{sp1}_{sp2}/vieref_{sp1}_{sp2}.{size}.step.windowed.weir.fst"
    params:
        output_format="raw_data/fst_step_window/pair_{sp1}_{sp2}/vieref_{sp1}_{sp2}.{size}.step",
        window_size="{size}",
    envmodules:
        "vcftools"
    shell:
        r"""
        vcftools --gzvcf {input.vcf} \
		--weir-fst-pop {input.sp1_sample} \
		--weir-fst-pop {input.sp2_sample} \
		--fst-window-size {params.window_size} \
		--out {params.output_format}
        """

rule r_plot_fst:
    input:
        rules.vcftools_step_window.output
    output:
        "results/fst_plots/pair_{sp1}_{sp2}/pair_{sp1}_{sp2}.{size}.fst.jpg"
    params:
        plot_title="Fst - {sp1} vs. {sp2}",
        window_info="window size: {size}bp"
    envmodules:
        "R"
    script:
        "../scripts/plot_fst.R"