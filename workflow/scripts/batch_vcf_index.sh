#!/bin/bash
#
#SBATCH --job-name=blastn
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=0-00:40:00
#SBATCH --partition=basic,short
#SBATCH --output=/lisc/scratch/botany/frschmidt/mk_test/logs/blastn_genes%j.log
#SBATCH --error=/lisc/scratch/botany/frschmidt/mk_test/logs/blastn_genes%j.err

ml bcftools

bcftools index -c --threads $SLURM_CPUS_PER_TASK data/variants/dio_nc_vieref.vcf.gz
