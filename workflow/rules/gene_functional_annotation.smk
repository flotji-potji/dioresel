rule subset_interpro:
    # rule subsets interpro results into 5 columns
    # it concatonates the GO terms
    # it sorts the accession description by length and takes longest name
    input:
        "data/interpro/transdecoder_pep_go_pathways.tsv"
    output:
        "raw_data/gene_annotation/interpro.tsv"
    shell:
        r"""
        awk -F'\t' -v OFS='\t' '{{print $1,$13,$14,$12,$4,length($13)}}' {input} \
            | sed -E 's/([0-9]*)\([A-Za-z]*\)/\1/g' \
            | sed 's/|/,/g' \
            | sed 's/\.t[0-9]//g' \
            | awk '$2 != "-" && $3 != "-"' \
            | sort -t$'\t' -k1,1 -k6,6 \
            | cut -f1,2,3,4,5 \
            | awk -F'\t' -v a='' -v lin='' \
            'a!=$1 {{
                if(NR!=1) {{ 
                    printf "%s\t", lin; 
                    for (g in gos) printf "%s,", g; 
                    printf "\t%s\t%s\t%s\n", $4, $5, "interpro"
                }}; 
                delete gos; split($3, gs,","); 
                for (i in gs) gos[gs[i]]++; a=$1; lin=$1"\t"$2
            }} 
            a==$1 && $3!="-" {{ 
                split($3,newg,","); 
                for (i in newg) if (!(newg[i] in gos)) gos[newg[i]]++
            }} 
            END {{
                printf "%s\t", lin; 
                for (g in gos) printf "%s,", g; 
                printf "\t%s\t%s\t%s\n", $4, $5, "interpro"
            }}' \
            | sed 's/\(.*\),/\1 /' \
            | sed 's/,-,/,/g' \
            > {output}
        """

rule subset_emapper:
    input:
        "data/emapper/emapper.annotations"
    output:
        "raw_data/gene_annotation/emapper.tsv"
    shell:
        r"""
        awk -F'\t' -v OFS='\t' '{{print $1,$8,$10,$2,$9,"eggnog"}}' {input} \
            | sed 's/\.t[0-9]//g' \
            | awk '$2 != "-" && $3 != "-"' \
            | sort -k1,1 | uniq \
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
        awk -F'\t' -v OFS='\t' \
            'FNR==NR{{a[$1]=$0}} 
            FNR!=NR {{ 
                if($1 in a) {{ 
                    if(a[$1] ~ "GO") {{
                        print a[$1]; 
                    }} else {{
                        print $0
                    }}
                    delete a[$1]
                }} else {{
                    print $0
                }}
            }} 
            END {{
                for (g in a) print a[g]
            }}' \
            <(sort -Vk1,1 {input.emapper}) \
            <(sort -Vk1,1 {input.interpro}) \
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
        
        awk -v OFS='\t' '{{print $0,"-","-","-","-","-","-"}}' $TMPDIR/tmp_out.gff >> {output}

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
            | sed -E 's/(ID=g[[:digit:]]*)([-g])/\1\t\2/g' \
            | cut -d$'\t' -f1,2,3,12,14,15,16,17,18 \
            | awk -F'\t' -v OFS='\t' '{{split($4, a, "="); $4=a[2]; print}}' \
            > {output}
        """
