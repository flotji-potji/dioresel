rule subset_interpro:
    input:
        "data/interpro/transdecoder_pep_go_pathways.tsv"
    output:
        "raw_data/gene_annotation/interpro.tsv"
    shell:
        r"""
        cut -f 1,13,14 {input} \
            | sed -E 's/([0-9]*)\([A-Za-z]*\)/\1/g' \
            | sed 's/|/,/g' \
            | sed 's/\.t[0-9]//g' \
            | awk '$2 != "-" && $3 != "-"' \
            > {output}
        """

rule subset_emapper:
    input:
        "data/emapper/emapper.annotations"
    output:
        "raw_data/gene_annotation/emapper.tsv"
    shell:
        r"""
        cut -f 1,8,10 {input} \
            | sed 's/\.t[0-9]//g' \
            | awk '$2 != "-" && $3 != "-"' \
            > {output}
        """

rule merge_emapper_interpro:
    input:
        emapper = rules.subset_emapper.output,
        interpro = rules.subset_interpro.output
    output:
        "raw_data/gene_annotation/func_gene_annotation.tsv"
    shell:
        r"""
        cat {input} \
            | sort -t$'\t' -Vk1,1 -k3,3r \
            | awk -F'\t' -v a='' 'a!=$1{{a = $1; print}}' \
            > {output}
        """

rule func_annotation_to_gff:
    input:
        func_anno = rules.merge_emapper_interpro.output,
        gff = "data/reference/vieillardii.gff"
    output:
        "raw_data/gene_annotation/func_gene_annotation.gff"
    shell:
        r"""
        awk 'FNR==NR{{a[$1]++}}; FNR!=NR{{split($9, b, "="); if(b[2] in a) print $0}}' \
            <(sort -k1 {input.func_anno}) \
            <(awk '$3~"gene"' {input.gff} | sort -k9) \
            > $TMPDIR/tmp_in.gff

        awk 'FNR==NR{{a[$1]++}}; FNR!=NR{{split($9, b, "="); if(!(b[2] in a)) print $0}}' \
            <(sort -k1 {input.func_anno}) \
            <(awk '$3~"gene"' {input.gff} | sort -k9) \
            > $TMPDIR/tmp_out.gff

        paste -d '\t' <(sort -k 9 -V $TMPDIR/tmp_in.gff) \
            <(sort -k1 -V {input.func_anno}) \
            > {output}
        
        awk -v OFS='\t' '{{print $0,"-","-","-"}}' $TMPDIR/tmp_out.gff >> {output}

        find $TMPDIR -maxdepth 1 -type f -delete
        """

rule func_annotation_to_window_bed:
    input:
        func_anno = rules.func_annotation_to_gff.output,
        windows = "data/reference/vieillardii.windows.bed"
    output:
        "results/gene_annotation/func_gene_annotation.windows.bed"
    shell:
        r"""
        bedtools intersect -a {input.windows} -b {input.func_anno} -wao \
            | grep -v '\-1' \
            | sed -E 's/([g][[:digit:]]*)([-g])/\1\t\2/g' \
            | cut -d$'\t' -f1,2,3,13,14,15 \
            > {output}
        """
