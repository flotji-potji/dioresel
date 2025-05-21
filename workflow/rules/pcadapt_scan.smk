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

rule bedtools_make_windows:
    input:
        "data/reference/vieillardii.fai"
    output:
        "data/reference/vieillardii.windows.bed"
    params:
        window_size = 10000
    shell:
        """
        bedtools makewindows -g {input} -w {params.window_size} > {output}
        """

rule pcadapt_to_bed:
    input:
        rules.pcadapt_scan_and_plot.output
    output:
        "raw_data/pcadapt/pair_{sp1}_{sp2}/single/pcadapt_snp_results.bed"
    shell:
        r"""
        awk 'BEGIN{{OFS="\t"}} NR>1{print $1, $2, $2+1, $5, $6}' {input} > {output}
        """

rule pcadapt_make_windows:
    input:
        windows = rules.bedtools_make_windows.output,
        pcadapt = rules.pcadapt_to_bed.output
    output:
        "raw_data/pcadapt/pair_{sp1}_{sp2}/single/pcadapt_snp_results.windowed.bed"
    shell:
        r"""
        pcadapt_length=$(wc -l {input.pcadapt} | cut -d' ' -f1)
        n_chunks=$(python -c "print(round($pcadapt_length/400000))")

        bedtools split -i {input.pcadapt} -n $n_chunks \
                       -p $TMPDIR/pcadapt_chunk -a simple

        for (file in $TMPDIR/pcadapt_chunk*); do
            bedtools intersect -a {input.windows} \
                               -b $file -wo > $TMPDIR/tmp.bed
        done

        bedtools groupby -g 1,2,3 -c 7,8,9 \
                         -o mean,mean,sum -i $TMPDIR/tmp.bed > $TMPDIR/group.bed 
        
        awk 'BEGIN{OFS="\t"} {print $1, $2+1, $3, $4, $5, $6}' $TMPDIR/group.bed > {output}

        rm $TMPDIR/*
        """

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