#!/bin/bash

# load modules
ml bcftools
ml conda
ml R
conda activate plink

# major part
# declare helper variables
VCF=$1
OUT=$(basename $1)
OUT=${OUT}_$(basename $2)
SCRIPTS=/lisc/project/botany/frschmidt/scripts

# subset vcf
bcftools view -Oz -S $2 $VCF -o $OUT.vcf.gz

# create bed file
plink --vcf $OUT.vcf.gz --double-id --allow-extra-chr \
		--make-bed --out $OUT
module unload conda

# do pcadapt plots
R --vanilla -s -f $SCRIPTS/plot_pcadapt.R --args $OUT $2
