rule annotate_vcf:
    input:
        vcf = "data/variants/pair_{sp1}_{sp2}/vieref_{sp1}_{sp2}.vcf.gz"
    output:
        annotated_vcf = "raw_data/snpeff_annotation/pair_{sp1}_{sp2}/vieref_{sp1}_{sp2}.smk.new.ann.vcf.gz"
    params:
        snpeff_prog = "/lisc/scratch/botany/frschmidt/bin/snpEff/snpEff.jar",
        snpeff_db = "/lisc/scratch/botany/frschmidt/mk_test/data/snpEff_DB/",
        snpeff_config = "/lisc/scratch/botany/frschmidt/bin/snpEff/snpEff.config",
        vieref_annot_version = "vie.refAnnot",
        out_dir = "raw_data/snpeff_annotation/pair_{sp1}_{sp2}"
    resources:
        mem_gb = 8,
        cpus_per_task = 2
    shell:
        """
        #java -Xmx8g -jar {params.snpeff_prog} {params.vieref_annot_version} \
        #        -dataDir {params.snpeff_db} \
        #        {input.vcf} | bgzip -c > {output.annotated_vcf}
        snpEff {params.vieref_annot_version} -c {params.snpeff_config} \
                -dataDir {params.snpeff_db} -noStats \
                {input.vcf} | bgzip -c > {output.annotated_vcf}

        """

