rule vcftools_sliding_window:
    # rule produces sliding window results as output
    input:
        vcf="data/variants/pair_{sp1}_{sp2}/vieref_{sp1}_{sp2}.vcf.gz",
        sp1_sample="data/variants/pair_{sp1}_{sp2}/{sp1}.samples",
        sp2_sample="data/variants/pair_{sp1}_{sp2}/{sp2}.samples"
    output:
        "raw_data/vcftools/fst_sliding_window/vieref_{sp1}_{sp2}.{size}.{step}.sliding.windowed.weir.fst"
    params:
        output_format="raw_data/vcftools/fst_sliding_window/vieref_{sp1}_{sp2}.{size}.{step}.sliding",
        window_size="{size}",
        window_step="{step}"
    envmodules:
        "vcftools"
    shell:
        r"""
        vcftools --gzvcf {input.vcf} \
		--weir-fst-pop {input.sp1_sample} \
		--weir-fst-pop {input.sp2_sample} \
		--fst-window-size {params.window_size} \
		--fst-window-step {params.window_step} \
		--out {params.output_format}
        """

rule r_plot_fst:
    input:
        rules.vcftools_sliding_window.output
    output:
        "results/fst_plots/pair_{sp1}_{sp2}.{size}.{step}.fst.jpg"
    params:
        plot_title="Fst - {sp1} vs. {sp2}",
        window_info="window size: {size}bp, window step: {step}bp"
    envmodules:
        "R"
    script:
        "../scripts/plot_fst.R"