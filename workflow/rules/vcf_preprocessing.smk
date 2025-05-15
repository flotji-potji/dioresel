rule bcftools_all_samples:
    input: 
        "data/variants/dio_nc_sub_vieref.vcf.gz"
    output:
        "data/variants/dio_nc_sub_vieref.all.samples"
    shell:
        """
        bcftools query -l {input} > {output}
        """

rule bcftools_subset_sample_pairs:
    input:
        rules.bcftools_all_samples.output
    output:
        sp1_sample = "data/variants/pair_{sp1}_{sp2}/{sp1}.samples",
        sp2_sample = "data/variants/pair_{sp1}_{sp2}/{sp2}.samples"
    params:
        sp1 = "{sp1}",
        sp2 = "{sp2}"
    shell:
        """
        grep {params.sp1} {input} > {output.sp1_sample}
        grep {params.sp2} {input} > {output.sp2_sample}
        """

rule bcftools_subset_vcf:
    input:
        vcf = rules.bcftools_all_samples.input,
        sp1_sample = rules.bcftools_subset_sample_pairs.output.sp1_sample,
        sp2_sample = rules.bcftools_subset_sample_pairs.output.sp2_sample
    output:
        "data/variants/pair_{sp1}_{sp2}/vieref_{sp1}_{sp2}.unfiltered.vcf.gz"
    shell:
        """
        bcftools view -S <(cat {input.sp1_sample} {input.sp2_sample}) \
                      -Oz -o {output} {input.vcf}
        """

rule vcftools_missing:
    input: 
        "data/variants/pair_{sp1}_{sp2}/vieref_{sp1}_{sp2}.unfiltered.vcf.gz"
    output:
        "data/variants/pair_{sp1}_{sp2}/vieref_{sp1}_{sp2}.missing.vcf.gz"
    shell:
        """
        vcftools --gzvcf {input} --max-missing 0.1 \
                 --recode --stdout | bgzip -c > {output}
        """

rule vcftools_maf:
    input:
        rules.vcftools_missing.output
    output:
        "data/variants/pair_{sp1}_{sp2}/vieref_{sp1}_{sp2}.maf.vcf.gz"
    shell:
        '''
        num_ind=$(bcftools query -l {input} | wc -l)
        maf=$(python -c "print(round(2/$num_ind, 3))")
        vcftools --gzvcf {input} --maf $maf \
                 --recode --stdout | bgzip -c > {output}
        '''

rule bcftools_monomorphic: 
    input:
        rules.vcftools_maf.output
    output:
        "data/variants/pair_{sp1}_{sp2}/vieref_{sp1}_{sp2}.vcf.gz"
    shell:
        """
        bcftools view -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' \
                      -Oz -o {output} {input}
        """