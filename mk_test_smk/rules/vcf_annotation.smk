rule snpeff_annotate_vcf:
    input:
        vcf = "data/variants/pair_{sp1}_{sp2}/vieref_{sp1}_{sp2}.vcf.gz"
    output:
        "raw_data/snpeff_annotation/pair_{sp1}_{sp2}/vieref_{sp1}_{sp2}.ann.vcf.gz"
    params:
        snpeff_db = "/lisc/scratch/botany/frschmidt/mk_test/data/snpEff_DB/",
        snpeff_config = "/lisc/scratch/botany/frschmidt/bin/snpEff/snpEff.config",
        vieref_annot_version = "vie.refAnnot",
        annot_stat = "raw_data/snpeff_annotation/pair_{sp1}_{sp2}/vie_annot_stat"
    shell:
        """
        snpEff {params.vieref_annot_version} -c {params.snpeff_config} \
                -dataDir {params.snpeff_db} -stats {params.annot_stat} \
                {input.vcf} | bgzip -c > {output}
        """

