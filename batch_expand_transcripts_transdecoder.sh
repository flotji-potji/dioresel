#!/bin/bash
#
#SBATCH --job-name=expTran
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=5G
#SBATCH --time=0-01:00:00
#SBATCH --partition=basic,short
#SBATCH --output=/lisc/scratch/botany/frschmidt/mk_test/logs/exp_trans_transdecoder%j.log
#SBATCH --error=/lisc/scratch/botany/frschmidt/mk_test/logs/exp_trans_ttransdecoder%j.err
#SBATCH --mail-type=END

# ENVIRONMENT 
ml transdecoder
ml bedtools

module list

# VARIABLES
WD=/lisc/scratch/botany/frschmidt/mk_test

# vieillardii reference fasta and annotation
REF_DIR=${WD}/data/reference
REF_FASTA=${REF_DIR}/vie_reference_genome/vieillardii1167c.asm.bp.p_ctg.fa
REF_FAI=${REF_DIR}/vie_reference_genome/vieillardii1167c.asm.bp.p_ctg.fa.fai
REF_ANNOT=${REF_DIR}/vieillardii_braker_tk.gff 

# output directories and files
raw_out_dir=${WD}/raw_data/reference
out_dir=${WD}/results/expanded_transcripts_transdecoder
gene_gff=${raw_out_dir}/vie_gene.gff # vie-gff gene file 
gene_gff_slop=${raw_out_dir}/vie_gene_slop1000bp.gff # vie-gff gene file with padded intervals
gene_fa=${raw_out_dir}/vie_gene_slop1000bp.fa # fasta-file of vie genes expanded 1000bp


# EXECUTION
mkdir -p ${raw_out_dir} && mkdir -p ${out_dir}

# First gene regions are extracted and intervals are padded to 1000bp
# Genes from the reference gff-file annotation get extracted
awk '$3 ~ /gene/ {print $0}' ${REF_ANNOT} > ${gene_gff}

# Genic regions are then expanded (padded) by 1000bp (-b) according to their strandedness (-s)
bedtools slop -i ${gene_gff} -g ${REF_FAI} -b 1000 -s > ${gene_gff_slop}

# extract padded genic regions to fasta file
bedtools getfasta -fi ${REF_FASTA} -bed ${gene_gff_slop} > ${gene_fa}

# Next TransDecoder predicts ORFs and CDS of gene transcripts
# LongOrfs extracts the longest ORFs from each genic region
TransDecoder.LongOrfs -t ${gene_fa} --output_dir ${out_dir}

# Predict decides based on the longest ORFs which to choose
TransDecoder.Predict -t ${gene_fa} --output_dir ${out_dir}


