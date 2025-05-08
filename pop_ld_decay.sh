#!/bin/bash

# load modules
ml conda
ml R
ml bcftools
conda activate plink

# major part
# declare helper variables
VCF=$1
SCRIPTS=/lisc/project/botany/frschmidt/scripts
SPP=$2

for file in $SPP/*; do
		OUT=$(basename $VCF)
		OUT=${OUT%.vcf.gz}	
		OUT=${OUT}_$(basename $file)

		# extract species from vcf
		bcftools view -Oz -S $file $VCF -o $OUT.vcf.gz

		# calculate ld decay with plink
		plink --vcf $OUT.vcf.gz --double-id --allow-extra-chr \
				--maf 0.01 --geno 0.1 --mind 0.5 \
				--thin 0.1 --r2 gz --ld-window 100 --ld-window-kb 1000 \
				--ld-window-r2 0 --out $OUT

		module unload conda

		# calculate average ld decay
		$SCRIPTS/ld_decay_calc.py -i $OUT.ld.gz -o $OUT

		# plot ld decay
		R --vanilla -s -f $SCRIPTS/plot_ld_decay.R --args $OUT.ld_decay_bins
done
