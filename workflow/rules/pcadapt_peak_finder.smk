rule manhattan_harvester_find_peaks:
    input:
        rules.pcadapt_scan_and_plot.output
    output:
        "raw_data/harvester/pair_{sp1}_{sp2}/found_peaks.out"
    params:
        harvester_prog="bin/harvester",
        out_dir="raw_data/harvester/pair_{sp1}_{sp2}",
        pval_column=5
    shell:
        r"""
        test -e harvester || cp {params.harvester_prog} ./
        mkdir -p {params.out_dir}
        ./harvester -chrcolumn 1 -lcolumn 2 -pcolumn {params.pval_column} \
                    -header yes -file {input} \
                    -out {output} -delim tab -missing NA \
                    -inlimit 0.001 -peak-limit 9 -dots 5 -shrink 2 || true
        test -e harvester && rm harvester
        """

rule manhattan_harvester_find_peaks_padj:
    input:
        rules.pcadapt_scan_and_plot.output
    output:
        "raw_data/harvester/pair_{sp1}_{sp2}/found_peaks_padj.out"
    params:
        harvester_prog="bin/harvester",
        out_dir="raw_data/harvester/pair_{sp1}_{sp2}",
        pval_column=6
    shell:
        r"""
        cp {params.harvester_prog} ./
        mkdir -p {params.out_dir}
        ./harvester -chrcolumn 1 -lcolumn 2 -pcolumn {params.pval_column} \
                    -header yes -file {input} \
                    -out {output} -delim tab -missing NA \
                    -inlimit 0.001 -peak-limit 5 -dots 5 -shrink 2 || true
        test -e harvester && rm harvester
        """

rule pcadapt_peak_finder:
    input:
        rules.manhattan_harvester_find_peaks.output
    output:
        "raw_data/harvester/pair_{sp1}_{sp2}/snps_of_interest.tsv"
    params:
        pcadapt_res=rules.pcadapt_scan_and_plot.output,
        raw_output_prefix="raw_data/harvester/pair_{sp1}_{sp2}/pval",
        plot_title="{sp1} - {sp2}",
        pval_threshold=0.05,
        pval_type="pval"
    script:
        "../scripts/pcadapt_peak_finder.R"

rule pcadapt_peak_finder_padj:
    input:
        rules.manhattan_harvester_find_peaks_padj.output
    output:
        "raw_data/harvester/pair_{sp1}_{sp2}/snps_of_interest_padj.tsv"
    params:
        pcadapt_res=rules.pcadapt_scan_and_plot.output,
        raw_output_prefix="raw_data/harvester/pair_{sp1}_{sp2}/padj",
        plot_title="{sp1} - {sp2}",
        pval_threshold=0.01,
        pval_type="padj"
    script:
        "../scripts/pcadapt_peak_finder.R"