rule bcftools_subset_cds:
    input:
        rules.snpeff_annotate_vcf.output
    output:
        "raw_data/snpeff_annotation/pair_{sp1}_{sp2}/vieref_{sp1}_{sp2}.ann.miss-syno.vcf.gz"
    shell:
        """
        bcftools view -i 'ANN ~ "missense" || ANN ~ "synonymous"' \
                      {input} -Oz -o {output} 
        """

rule vcftools_pop_freq:
    input:
        vcf = rules.bcftools_subset_cds.output,
        sample1 = "data/variants/pair_{sp1}_{sp2}/{sp1}.samples",
        sample2 = "data/variants/pair_{sp1}_{sp2}/{sp2}.samples"
    output:
        sample1_out = "raw_data/pop_allel_freq/pair_{sp1}_{sp2}/vieref_{sp1}_{sp2}.sp_{sp1}.frq",
        sample2_out = "raw_data/pop_allel_freq/pair_{sp1}_{sp2}/vieref_{sp1}_{sp2}.sp_{sp2}.frq"
    params:
        out_pre_sample1 = "raw_data/pop_allel_freq/pair_{sp1}_{sp2}/vieref_{sp1}_{sp2}.sp_{sp1}",
        out_pre_sample2 = "raw_data/pop_allel_freq/pair_{sp1}_{sp2}/vieref_{sp1}_{sp2}.sp_{sp2}" 
    shell:
        """
        vcftools --gzvcf {input.vcf} --keep {input.sample1} \
                 --freq --out {params.out_pre_sample1} 
        vcftools --gzvcf {input.vcf} --keep {input.sample2} \
                 --freq --out {params.out_pre_sample2} 
        """