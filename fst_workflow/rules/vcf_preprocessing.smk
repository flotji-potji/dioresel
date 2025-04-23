rule vcftools_missing:
    input: 
        "data/variants/pair_{sp1}_{sp2}/vieref_{sp1}_{sp2}.unfiltered.vcf.gz"
    output:
        "data/variants/pair_{sp1}_{sp2}/vieref_{sp1}_{sp2}.missing.vcf.gz"
    resources:
        mem_gb = 8,
        cpus_per_task = 2
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
    resources:
        mem_gb = 8,
        cpus_per_task = 2
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
    resources:
        mem_gb = 8,
        cpus_per_task = 2
    shell:
        """
        bcftools view -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' \
                      -Oz -o {output} {input}
        """