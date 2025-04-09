#!/bin/bash
#
#SBATCH --job-name=transdecS
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=5G
#SBATCH --time=0-01:00:00
#SBATCH --partition=basic,short
#SBATCH --output=/lisc/scratch/botany/frschmidt/mk_test/logs/transdecoder_singleBest%j.log
#SBATCH --error=/lisc/scratch/botany/frschmidt/mk_test/logs/transdecoder_singleBest%j.err
#SBATCH --mail-type=END

# ENVIRONMENT 
ml transdecoder

module list

# VARIABLES

WD=/lisc/scratch/botany/frschmidt/mk_test

# output directories and files
raw_out_dir=${WD}/raw_data/reference
out_dir=${WD}/raw_data/transdecoder_single_best
gene_fa=${raw_out_dir}/vie_exon.fa

# EXECUTION
mkdir -p ${raw_out_dir} && mkdir -p ${out_dir}

# Next TransDecoder prediction
TransDecoder.LongOrfs --complete_orfs_only -t ${gene_fa} --output_dir ${out_dir}

TransDecoder.Predict -t ${gene_fa} --single_best_only --output_dir ${out_dir} 

