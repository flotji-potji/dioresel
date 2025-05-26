rule filter_fst:
    input:
        rules.vcftools_to_bed.output
    output:
        "raw_data/filtered_bed/pair_{sp1}_{sp2}/fst.filtered.bed"
    shell:
        r"""
        awk '($5 > 3 && $4 > 0){{print}}' {input} > {output}
        """

rule outliers_fst:
    input:
        rules.filter_fst.output
    output:
        "raw_data/filtered_bed/pair_{sp1}_{sp2}/fst.outliers.bed"
    params:
        prob = 0.95,
        col = 4,
        sign = ">"
    shell:
        r"""
        quantile=$(cut -f {params.col} {input} \
            | Rscript -e 'print(quantile(scan("stdin", quiet=TRUE), {params.prob}))' \
            | tail -n +2)
        awk -v quant=$quantile '${params.col} {params.sign} quant {{print}}' {input} > {output}
        """

rule filter_pi:
    input:
        rules.pi_to_bed.output
    output:
        "raw_data/filtered_bed/pair_{sp1}_{sp2}/pi.filtered.bed"
    shell:
        r"""
        awk '($5 > 3 && $4 != "nan"){{print}}' {input} > {output}
        """

use rule filter_pi as filter_pair_pi with:
    input:
        rules.pair_pi_to_bed.output
    output:
        "raw_data/filtered_bed/pair_{sp1}_{sp2}/pi_pair.filtered.bed"

use rule outliers_fst as outliers_pi with:
    input:
        rules.filter_pi.output
    output:
        "raw_data/filtered_bed/pair_{sp1}_{sp2}/pi.outliers.bed"

use rule outliers_fst as outliers_pi_bottom with:
    input:
        rules.filter_pi.output
    output:
        "raw_data/filtered_bed/pair_{sp1}_{sp2}/pi.outliers.bottom.bed"
    params:
        prob = 0.05,
        sign = "<"

use rule outliers_pi as outliers_pair_pi with:
    input:
        rules.filter_pair_pi.output
    output:
        "raw_data/filtered_bed/pair_{sp1}_{sp2}/pi_pair.outliers.bed"

use rule outliers_pi_bottom as outliers_pair_pi_bottom with:
    input:
        rules.filter_pair_pi.output
    output:
        "raw_data/filtered_bed/pair_{sp1}_{sp2}/pi_pair.outliers.bottom.bed"

use rule filter_pi as filter_tajimad with:
    input:
        rules.tajimad_to_bed.output
    output:
        "raw_data/filtered_bed/pair_{sp1}_{sp2}/tajimad.filtered.bed"

use rule filter_tajimad as filter_pair_tajimad with:
    input:
        rules.pair_tajimad_to_bed.output
    output:
        "raw_data/filtered_bed/pair_{sp1}_{sp2}/tajimad_pair.filtered.bed"

use rule outliers_fst as outliers_sel_tajimad with:
    input:
        rules.filter_tajimad.output
    output:
        "raw_data/filtered_bed/pair_{sp1}_{sp2}/tajimad.sel_outliers.bed"
    params:
        prob = 0.05,
        sign = "<"
        
use rule outliers_sel_tajimad as outliers_sel_pair_tajimad with:
    input:
        rules.filter_pair_tajimad.output
    output:
        "raw_data/filtered_bed/pair_{sp1}_{sp2}/tajimad_pair.sel_outliers.bed"

rule filter_pcadapt:
    input:
        rules.pcadapt_make_windows.output
    output:
        "raw_data/filtered_bed/pair_{sp1}_{sp2}/pcadapt.filtered.bed"
    shell:
        r"""
        awk '($6 > 3 && $5 != "NA" && $4 != "NA"){{print}}' {input} > {output}
        """

rule outliers_pcadapt:
    input:
        rules.filter_pcadapt.output
    output:
        "raw_data/filtered_bed/pair_{sp1}_{sp2}/pcadapt.outliers.bed"
    params:
        prob = 0.05
    shell:
        r"""
        awk '$4 < {params.prob} {{print}}' {input} > {output}
        """

rule outliers_sel_mk_test:
    input:
        rules.mk_test_to_window_bed.output
    output:
        "raw_data/filtered_bed/pair_{sp1}_{sp2}/mk_test.outliers.bed"
    params:
        prob = 0.05
    shell:
        r"""
        awk '($6 < {params.prob} && $5 > $4) {{print}}' {input} > {output}
        """
