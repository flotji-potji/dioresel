rule plink_convert_bed:
    input:
        "data/variants/pair_{sp1}_{sp2}/vieref_{sp1}_{sp2}.vcf.gz"
    output:
        "data/variants/pair_{sp1}_{sp2}/vieref_{sp1}_{sp2}.bed"
    params:
        output_prefix="data/variants/pair_{sp1}_{sp2}/vieref_{sp1}_{sp2}"
    shell:
        r"""
        plink --vcf {input} --double-id --allow-extra-chr \
              --make-bed --out {params.output_prefix}
        """

rule pcadapt_scan_and_plot:
    input:
        rules.plink_convert_bed.output
    output:
        "raw_data/pcadapt/pair_{sp1}_{sp2}/single/pcadapt_snp_results.tsv"
    params:
        input_prefix=rules.plink_convert_bed.params.output_prefix,
        raw_output_prefix="raw_data/pcadapt/pair_{sp1}_{sp2}/single/",
        k_param=1,
        plot_title="{sp1} - {sp2}",
        ld_thin="no",
        vie_genome_size=1523932088
    resources:
        mem_mb=4000,
        cpus_per_task=2
    script:
        "../scripts/pcadapt_scan_and_plot.R"

rule pcadapt_ldthin_scan_and_plot:
    input:
        rules.plink_convert_bed.output
    output:
        "raw_data/pcadapt/pair_{sp1}_{sp2}/ldthin/pcadapt_snp_results.tsv"
    params:
        input_prefix=rules.plink_convert_bed.params.output_prefix,
        raw_output_prefix="raw_data/pcadapt/pair_{sp1}_{sp2}/ldthin/",
        k_param=1,
        plot_title="LD-thin {sp1} - {sp2}",
        ld_thin="yes",
        vie_genome_siz=1523932088
    resources:
        mem_mb=4000,
        cpus_per_task=2
    script:
        "../scripts/pcadapt_scan_and_plot.R"