#!/bin/bash
#
#SBATCH --job-name=pixy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=8G
#SBATCH --time=0-00:40:00
#SBATCH --partition=basic,short
#SBATCH --output=/lisc/scratch/botany/frschmidt/dioresel/logs/pixy%j.log
#SBATCH --error=/lisc/scratch/botany/frschmidt/dioresel/logs/pixy%j.err

ml conda
conda activate pixy

files=(data/pixy_pops/*.pop)

file=${files[$SLURM_ARRAY_TASK_ID]}

pixy --stats dxy \
	 --vcf data/variants/dio_nc_vieref.vcf.gz \
	 --populations $file \
	 --window_size 10000  \
	 --bypass_invariant_check \
	 --output_folder raw_data/pixy/ \
	 --output_prefix $(basename $file) \
	 --n_cores $SLURM_CPUS_PER_TASK

	 #--chromosomes 'ptg000001l' \


