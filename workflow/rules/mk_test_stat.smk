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

rule mk_test_to_gff:
    input:
        mk_test = rules.mk_test.output[1],
        gff = "data/reference/vieillardii.gff"
    output:
        "raw_data/mk_test/pair_{sp1}_{sp2}/mkt_table.gff"
    shell:
        r"""
        awk 'FNR==NR{{a[$1]++}}; FNR!=NR{{split($9, b, "="); if(b[2] in a) print $0}}' \
            <(tail -n +2 {input.mk_test} | sort -k1) \
            <(awk '$3~"gene"' {input.gff} | sort -k9) \
            > $TMPDIR/tmp.gff

        paste $TMPDIR/tmp.gff <(tail -n +2 {input.mk_test} | sort -k1) \
            > {output}

        find $TMPDIR -maxdepth 1 -type f -delete
        """


rule mk_test_to_window_bed:
    input:
        mk_test = rules.mk_test_to_gff.output,
        windows = "data/reference/vieillardii.windows.bed"
    output:
        "raw_data/mk_test/pair_{sp1}_{sp2}/mkt_table.window.bed"
    shell:
        r"""
        bedtools intersect -a {input.windows} -b {input.mk_test} -wao \
            | grep -v '\-1' \
            | cut -f 1,2,3,17,18,19,20 \
            > {output}
        """

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