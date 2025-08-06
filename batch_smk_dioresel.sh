#!/bin/bash
#
#SBATCH --job-name=smk_dioresel
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=0-13:00:00
#SBATCH --output=/lisc/scratch/botany/frschmidt/dioresel/logs/smk_dioresel%j.log
#SBATCH --error=/lisc/scratch/botany/frschmidt/dioresel/logs/smk_dioresel%j.err

### ENVIRONMENT
ml conda
conda activate dioresel

module list

### EXECUTION
snakemake -c 1 --workflow-profile config/slurm --use-conda 
