import os.path

rule merge_exons:
    input:
        ref_fasta="data/reference/vie_reference_genome/vieillardii1167c.asm.bp.p_ctg.fa",
        ref_gff="data/reference/vieillardii_braker_tk.gff"
    output:
        exon_fasta="raw_data/reference/vie_exon.fa"
    params:
        out_dir="raw_data/reference",
        exon_fa="vie_exon.fa",
        tmp_fasta="raw_data/reference/vie.fa"
    conda:
        "agat"
    shell:
        """
        fold {input.ref_fasta} > {params.tmp_fasta}
        cd {params.out_dir}
        agat_sp_extract_sequences.pl --gff {input.ref_gff} --fasta {params.tmp_fasta} -t exon --merge -o {output.exon_fasta}
        """

rule longest_orf:
    input:
        transcripts=rules.merge_exons.output.exon_fasta
    output:
        transdecoder_dir="raw_data/transdecoder"
    envmodules:
        "transdecoder"
    shell:
        """
        TransDecoder.LongOrfs --complete_orfs_only -t {input.transcripts} --output_dir {output.transdecoder_dir}
        """

rule pep_chunks:
    input:
        pep_fasta="{rules.longest_orf.output}/{os.path.basename(rules.merge_exons.output)}.transdecoder_dir/longest_orf.pep"
    output:
        expand("raw_data/blast_chunks/chunk{n}.pep.fa", n=range(1, 26))    
    params:
        tmp_fa="{input.pep_fasta}.fa",
        out_chunk="raw_data/blast_chunks/chunk"
    envmodules:
        "samtools"
    shell:
        r"""
        sed 's/\*//g' {input} > {params.tmp_fa}
        samtools faidx {tmp_fa}
        list=($(cut -f1 {tmp_fa}.fai))
        for i in {0..26}; do
            si=$(($i*1000))
            samtools faidx {tmp_fa} ${list[@]:$si:1000} > {params.out_chunk}$(($i+1)).pep.fa
        done
        """

rule blast_orfs:
    input:
        query="{rules.pep_chunks.params.out_chunk}{n}.pep.fa"
    output:
        "raw_data/blastp/chunk{n}.blastp.outfmt6"
    params:
        db="/lisc/scratch/mirror/ncbi/2025-02-25/nr"
    envmodules:
        "ncbiblastplus"
    shell:
        """
        blastp -query {input.query} -db {params.db} -out {output.blast_query} -num_threads {$SLURM_CPUS_PER_TASK} -outfmt 6 -max_target_seqs 1 -evalue 1e-5
        """