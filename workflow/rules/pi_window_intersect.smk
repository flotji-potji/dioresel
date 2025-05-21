rule filter_fst:
    input:
        rules.vcftools_to_bed.output
    output:
        "raw_data/pi_intersect/pair_{sp1}_{sp2}/fst.filtered.bed"
    shell:
        r"""
        awk '($5 > 3 & $4 != "nan"){{print}}' {input} > {output}
        """

use rule filter_fst as filter_pi with:
    input:
        rules.pi_to_bed.output
    output:
        "raw_data/pi_intersect/pair_{sp1}_{sp2}/pi.filtered.bed"

use rule filter_fst as filter_tajimad with:
    input:
        rules.tajimad_to_bed.output
    output:
        "raw_data/pi_intersect/pair_{sp1}_{sp2}/tajimad.filtered.bed"

rule filter_pcadapt:
    input:
        rules.pcadapt_make_windows.output
    output:
        "raw_data/pi_intersect/pair_{sp1}_{sp2}/pcadapt.filtered.bed"
    shell:
        r"""
        awk '($6 > 3 & $5 != "NA" & $4 != "NA"){{print}}' {input} > {output}
        """

rule pi_top_windows:
    input:
        rules.filter_pi.output
    output:
        "raw_data/pi_intersect/pair_{sp1}_{sp2}/pi.top_windows.bed"
    params:
        num_windows = 5000,
        order = "-r"
    shell:
        r"""
        sort -g -k4,4 {params.order} {input} | head -n {params.num_windows} > {output}
        """

rule pi_window_intersect:
    input:
        pi = rules.pi_top_windows.output,
        tajimad = rules.filter_tajimad.output,
        pcadapt = rules.filter_pcadapt.output,
        fst = rules.filter_fst.output
    output:
        "raw_data/pi_intersect/pair_{sp1}_{sp2}/pi.tajimad.pcadapt.fst.bed"
    shell:
        r"""
        bedtools intersect -a {input.pi} \
                           -b {input.tajimad} {input.pcadapt} {input.fst} \
                           -wo > {output}
        """