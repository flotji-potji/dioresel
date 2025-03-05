#!/bin/bash
#
#SBATCH --job-name=tranD
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=5G
#SBATCH --time=0-01:00:00
#SBATCH --partition=basic,short
#SBATCH --output=/lisc/scratch/botany/frschmidt/mk_test/logs/transdecoder%j.log
#SBATCH --error=/lisc/scratch/botany/frschmidt/mk_test/logs/transdecoder%j.err
#SBATCH --mail-type=END

# ENVIRONMENT 
ml transdecoder
ml samtools

module list

# VARIABLES

WD=/lisc/scratch/botany/frschmidt/mk_test

# vieillardii reference fasta and annotation
REF_DIR=${WD}/data/reference
REF_FASTA=${REF_DIR}/vie_reference_genome/vieillardii1167c.asm.bp.p_ctg.fa
REF_ANNOT=${REF_DIR}/vieillardii_braker_tk.gff 

# output directories and files
raw_out_dir=${WD}/raw_data/reference
out_dir=${WD}/results/transdecoder
gene_reg=${raw_out_dir}/vie_gene.regions # contains region identifiers of genes
gene_fa=${raw_out_dir}/vie_gene.fa # fasta-file of vie genes


# EXECUTION
mkdir -p ${raw_out_dir} && mkdir -p ${out_dir}

# First the gene positions, in faidx format,
# need to be extracted from the gff-file
awk '{if ($3 == "gene") print $1":"$4"-"$5}' ${REF_ANNOT} > ${gene_reg}

# Next genes from reference assembly can be
# extracted with gene regions
samtools faidx ${REF_FASTA} -r ${gene_reg} -o ${gene_fa}

# Next TransDecoder prediction
TransDecoder.LongOrfs -t ${gene_fa} --output_dir ${out_dir}

TransDecoder.Predict -t ${gene_fa} --output_dir ${out_dir}


