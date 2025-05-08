#!/bin/bash
#SBATCH --job-name=test_pca
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --mail-type=ALL
#SBATCH --time=0-0:01:00
#SBATCH --output=log/%x-%j_%a.out
#SBATCH --error=log/%x-%j_%a.err
# sbatch time more like 10 sec

# load modules
ml conda
ml R
conda activate plink

# major part
# declare helper variables
VCF=$1
OUT=$(basename $1)
OUT=${OUT%.vcf.gz}
SCRIPTS=/lisc/project/botany/frschmidt/scripts

# calculate PCAs with plink
plink --vcf $VCF --double-id --allow-extra-chr \
	--pca --out $OUT

# plot PCAs
R --vanilla \
	-s -f $SCRIPTS/plot_pca.R \
	--args $OUT.eigenvec $OUT.eigenval $2 $3
