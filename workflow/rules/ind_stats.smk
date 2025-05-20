rule vcftools_pi_stat:
    input:
        vcf = "data/variants/pair_{sp1}_{sp2}/vieref_{sp1}_{sp2}.vcf.gz",
        sample = "data/variants/pair_{sp1}_{sp2}/{sp2}.samples"
    output:
        pi = "raw_data/ind_stat/pair_{sp1}_{sp2}/sp_{sp2}.windowed.pi"
    params:
        output_format = "raw_data/ind_stat/pair_{sp1}_{sp2}/sp_{sp2}",
        window_size = 10000
    shell:
        r"""
        vcftools --gzvcf {input.vcf} \
        --keep {input.sample} \
		--window-pi {params.window_size} \
		--out {params.output_format}
        """

use rule vcftools_to_bed as pi_to_bed with:
    input:
        rules.vcftools_pi_stat.output.pi
    output:
        "raw_data/ind_stat/pair_{sp1}_{sp2}/sp_{sp2}.windowed.pi.bed"
    params:
        header = "PI"

rule vcftools_pi_stat_filtered:
    input:
        vcf = "data/variants/pair_{sp1}_{sp2}/vieref_{sp1}_{sp2}.vcf.gz",
        sample = "data/variants/pair_{sp1}_{sp2}/{sp2}.samples"
    output:
        pi = "raw_data/ind_stat/pair_{sp1}_{sp2}/sp_{sp2}.filtered.windowed.pi"
    params:
        output_format = "raw_data/ind_stat/pair_{sp1}_{sp2}/sp_{sp2}.filtered",
        window_size = 10000
    shell:
        r"""
        num_ind=$(bcftools query -l {input} | wc -l)
        maf=$(python -c "print(round(2/$num_ind, 3))")

        vcftools --gzvcf {input.vcf} \
        --keep {input.sample} \
        --max-missing 0.1 \
        --maf $maf \
		--window-pi {params.window_size} \
		--out {params.output_format}
        """

use rule vcftools_to_bed as pi_filtered_to_bed with:
    input:
        rules.vcftools_pi_stat_filtered.output.pi
    output:
        "raw_data/ind_stat/pair_{sp1}_{sp2}/sp_{sp2}.filtered.windowed.pi.bed"
    params:
        header = "PI"

rule r_plot_pi:
    input:
        rules.vcftools_pi_stat.output.pi,
        rules.pi_to_bed.output
    output:
        "results/pi_plots/pair_{sp1}_{sp2}/sp_{sp2}.jpg"
    params:  
        plot_title = "PI - {sp1} vs. {sp2}",
        window_info = "window size: 10000bp",
        col_type = "PI"
    script:
        "../scripts/plot_vcftools_results.R"

rule r_plot_pi_filtered:
    input:
        rules.vcftools_pi_stat_filtered.output.pi,
        rules.pi_filtered_to_bed.output
    output:
        "results/pi_plots/pair_{sp1}_{sp2}/sp_{sp2}.filtered.jpg"
    params:  
        plot_title = "PI - {sp1} vs. {sp2}",
        window_info = "filtered; window size: 10000bp",
        col_type = "PI"
    script:
        "../scripts/plot_vcftools_results.R"

rule vcftools_tajima_stat:
    input:
        vcf = "data/variants/pair_{sp1}_{sp2}/vieref_{sp1}_{sp2}.vcf.gz",
        sample = "data/variants/pair_{sp1}_{sp2}/{sp2}.samples"
    output:
        tajimad = "raw_data/ind_stat/pair_{sp1}_{sp2}/sp_{sp2}.Tajima.D"
    params:
        output_format = "raw_data/ind_stat/pair_{sp1}_{sp2}/sp_{sp2}",
        window_size = 10000
    shell:
        r"""
        vcftools --gzvcf {input.vcf} \
        --keep {input.sample} \
        --TajimaD {params.window_size} \
		--out {params.output_format}
        """

use rule vcftools_to_bed as tajimad_to_bed with:
    input:
        rules.vcftools_tajima_stat.tajimad
    output:
        "raw_data/ind_stat/pair_{sp1}_{sp2}/sp_{sp2}.Tajima.D.bed"
    params:
        header = "TajimaD"

rule vcftools_tajima_stat_filtered:
    input:
        vcf = "data/variants/pair_{sp1}_{sp2}/vieref_{sp1}_{sp2}.vcf.gz",
        sample = "data/variants/pair_{sp1}_{sp2}/{sp2}.samples"
    output:
        tajimad = "raw_data/ind_stat/pair_{sp1}_{sp2}/sp_{sp2}.filtered.Tajima.D"
    params:
        output_format = "raw_data/ind_stat/pair_{sp1}_{sp2}/sp_{sp2}.filtered",
        window_size = 10000
    shell:
        r"""
        num_ind=$(bcftools query -l {input} | wc -l)
        maf=$(python -c "print(round(2/$num_ind, 3))")

        vcftools --gzvcf {input.vcf} \
        --keep {input.sample} \
        --max-missing 0.1 \
        --maf $maf \
        --TajimaD {params.window_size} \
		--out {params.output_format}
        """

use rule vcftools_to_bed as tajimad_filtered_to_bed with:
    input:
        rules.vcftools_tajima_stat.filtered.tajimad
    output:
        "raw_data/ind_stat/pair_{sp1}_{sp2}/sp_{sp2}.filtered.Tajima.D.bed"
    params:
        header = "TajimaD"

rule r_plot_tajimad:
    input:
        rules.vcftools_tajima_stat.output.tajimad,
        rules.tajimad_to_bed.output
    output:
        "results/tajimad_plots/pair_{sp1}_{sp2}/sp_{sp2}.jpg"
    params:  
        plot_title = "Tajima's D - {sp1} vs. {sp2}",
        window_info = "window size: 10000bp",
        col_type = "TajimaD"
    script:
        "../scripts/plot_vcftools_results.R"

rule r_plot_tajimad_filtered:
    input:
        rules.vcftools_tajima_stat_filtered.output.tajimad,
        rules.tajimad_filtered_to_bed.output
    output:
        "results/tajimad_plots/pair_{sp1}_{sp2}/sp_{sp2}.filtered.jpg"
    params:  
        plot_title = "Tajima's D - {sp1} vs. {sp2}",
        window_info = "filtered; window size: 10000bp",
        col_type = "TajimaD"
    script:
        "../scripts/plot_vcftools_results.R"