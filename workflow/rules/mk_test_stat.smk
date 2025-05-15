rule mkt_summary:
    input:
        vcf = rules.bcftools_subset_cds.output,
        sample1 = rules.vcftools_pop_freq.output.sample1_out, 
        sample2 = rules.vcftools_pop_freq.output.sample2_out
    output:
        gene_summary = "raw_data/mkt_summary/pair_{sp1}_{sp2}/gene_summary.Rdata",
        vcf_mkt = "raw_data/mkt_summary/pair_{sp1}_{sp2}/vcf_mkt_df.Rdata"
    params:
        output_prefix = "raw_data/mkt_summary/pair_{sp1}_{sp2}"
    script:
        "../scripts/mkt_summary.R"

rule mk_test:
    input:
        rules.mkt_summary.output.gene_summary
    output:
        "raw_data/mk_test/pair_{sp1}_{sp2}/mkt_table.Rdata",
        "raw_data/mk_test/pair_{sp1}_{sp2}/mkt_table.tsv"
    params:
        output_prefix = "raw_data/mk_test/pair_{sp1}_{sp2}"
    script:
        "../scripts/mk_test.R"

rule dnds_test:
    input:
        rules.mkt_summary.output.gene_summary
    output:
        "raw_data/dnds_test/pair_{sp1}_{sp2}/dnds_table.Rdata",
        "raw_data/dnds_test/pair_{sp1}_{sp2}/dnds_table.tsv"
    params:
        output_prefix = "raw_data/dnds_test/pair_{sp1}_{sp2}"
    script:
        "../scripts/dnds_test.R"