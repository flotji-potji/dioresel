#!/bin/bash
#
#SBATCH --job-name=smk_ind
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=0-03:00:00
#SBATCH --output=/lisc/scratch/botany/frschmidt/ind_stat/logs/smk_ind_stat%j.log
#SBATCH --error=/lisc/scratch/botany/frschmidt/ind_stat/logs/smk_ind_stat%j.err

### ENVIRONMENT
ml conda
conda activate diosel

module list

### VARIABLES
WD=/lisc/scratch/botany/frschmidt/ind_stat

# vieillardii reference fasta and annotation
smk_dir=${WD}/scripts/ind_stat_workflow
smk_profile=${smk_dir}/profiles
smk_workflow=${smk_dir}/Snakefile

#wildcard="raw_data/snpeff_annotation/pair_fla_spn/vieref_fla_spn.smk.ann.vcf.gz"

### EXECUTION
snakemake --cores 1 --profile ${smk_profile} \
          -s ${smk_workflow} --use-conda ${wildcard}
