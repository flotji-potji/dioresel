rule vcftools_pi_pair_stat:
    input:
        vcf = "data/variants/pair_{sp1}_{sp2}/vieref_{sp1}_{sp2}.vcf.gz",
    output:
        pi = "raw_data/pair_stat/pair_{sp1}_{sp2}/pair_{sp1}_{sp2}.windowed.pi"
    params:
        output_format = "raw_data/pair_stat/pair_{sp1}_{sp2}/pair_{sp1}_{sp2}",
        window_size = 10000
    shell:
        r"""
        vcftools --gzvcf {input.vcf} \
		--window-pi {params.window_size} \
		--out {params.output_format}
        """

rule r_plot_pair_pi:
    input:
        rules.vcftools_pi_pair_stat.output.pi
    output:
        "results/pi_plots/pair_{sp1}_{sp2}/pair_{sp1}_{sp2}.jpg"
    params:  
        plot_title = "PI - {sp1} vs. {sp2}",
        window_info = "window size: 10000bp",
        col_type = "PI"
    script:
        "../scripts/plot_vcftools_results.R"

rule vcftools_tajima_pair_stat:
    input:
        vcf = "data/variants/pair_{sp1}_{sp2}/vieref_{sp1}_{sp2}.vcf.gz",
    output:
        tajimad = "raw_data/pair_stat/pair_{sp1}_{sp2}/pair_{sp1}_{sp2}.Tajima.D",
    params:
        output_format = "raw_data/pair_stat/pair_{sp1}_{sp2}/pair_{sp1}_{sp2}",
        window_size = 10000
    shell:
        r"""
        vcftools --gzvcf {input.vcf} \
        --TajimaD {params.window_size} \
		--out {params.output_format}
        """

rule r_plot_pair_tajimad:
    input:
        rules.vcftools_tajima_pair_stat.output.tajimad
    output:
        "results/tajimad_plots/pair_{sp1}_{sp2}/pair_{sp1}_{sp2}.jpg"
    params:  
        plot_title = "Tajima's D - {sp1} vs. {sp2}",
        window_info = "window size: 10000bp",
        col_type = "TajimaD"
    script:
        "../scripts/plot_vcftools_results.R"