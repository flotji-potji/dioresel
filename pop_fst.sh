#!/bin/bash
#SBATCH --job-name=pc_fst
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --mail-type=ALL
#SBATCH --time=0-0:02:00
#SBATCH --output=log/%x-%j_%a.out
#SBATCH --error=log/%x-%j_%a.err

# load modules
ml vcftools

# define standard variables
SCRIPTS=/lisc/project/botany/frschmidt/scripts
CHUNK=/lisc/scratch/botany/frschmidt/data/all_chr_cs/chunk${SLURM_ARRAY_TASK_ID}_*
CHUNK_SHORT=$(basename $CHUNK)
CHUNK_SHORT=${CHUNK_SHORT%.vcf.gz}
CLASS_DIR=/lisc/scratch/botany/frschmidt/data/individuals_cs
SPP_DIR=$CLASS_DIR/spp_class

# calculate fst with vcftools
cd $CHUNK_SHORT
mkdir -p fst
cd fst

vcftools --gzvcf $CHUNK \
		--weir-fst-pop $SPP_DIR/spicnga \
		--weir-fst-pop $SPP_DIR/calciphila \
		--out ./$CHUNK_SHORT
