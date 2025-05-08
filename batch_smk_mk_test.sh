#!/bin/bash
#
#SBATCH --job-name=smk-mk
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=2G
#SBATCH --time=0-03:00:00
#SBATCH --partition=basic,short
#SBATCH --output=/lisc/scratch/botany/frschmidt/mk_test/logs/smk_mktest%j.log
#SBATCH --error=/lisc/scratch/botany/frschmidt/mk_test/logs/smk_mktest%j.err

### ENVIRONMENT
ml conda
conda activate diosel

module list

### VARIABLES
WD=/lisc/scratch/botany/frschmidt/mk_test

# vieillardii reference fasta and annotation
smk_dir=${WD}/scripts/mk_test_smk
smk_profile=${smk_dir}/profiles
smk_workflow=${smk_dir}/rules/vcf_annotation.smk

wildcard="raw_data/snpeff_annotation/pair_fla_spn/vieref_fla_spn.smk.new.ann.vcf.gz"

### EXECUTION
snakemake --cores 1 --profile ${smk_profile} \
          -s ${smk_workflow} --use-conda ${wildcard}
