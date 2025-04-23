#!/bin/bash
#
#SBATCH --job-name=smk_fst
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --mail-type=END
#SBATCH --time=0-03:00:00
#SBATCH --output=/lisc/scratch/botany/frschmidt/fst/logs/smk_fst%j.log
#SBATCH --error=/lisc/scratch/botany/frschmidt/fst/logs/smk_fst%j.err

### ENVIRONMENT
ml conda
conda activate diosel

module list

### VARIABLES
WD=/lisc/scratch/botany/frschmidt/fst

# vieillardii reference fasta and annotation
smk_dir=${WD}/scripts/fst_workflow
smk_profile=${smk_dir}/profiles
smk_workflow=${smk_dir}/Snakefile

#wildcard="raw_data/snpeff_annotation/pair_fla_spn/vieref_fla_spn.smk.ann.vcf.gz"

### EXECUTION
snakemake --cores 1 --profile ${smk_profile} \
          -s ${smk_workflow} --use-conda ${wildcard}
