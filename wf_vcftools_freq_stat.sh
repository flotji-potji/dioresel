#!/bin/bash

# ENVIRONMENT 
ml vcftools

module list

# VARIABLES

WD=/lisc/scratch/botany/frschmidt/mk_test

# directories of stored files and output directories
raw_data_dir=${WD}/raw_data
pair_var_dir=${WD}/data/variants/pair_cal_spn/
pair_ann_dir=${raw_data_dir}/snpeff_annotation/pair_cal_spn
out_dir=${raw_data_dir}/pop_allel_freq/pair_cal_spn

# files used for execution
ann_vcf=${pair_ann_dir}/vieref_cal_spn.ann.miss-syno.vcf.gz
sp1_sample=${pair_var_dir}/calciphila.samples
sp2_sample=${pair_var_dir}/sppicnga.samples
sp1_out=${out_dir}/$(basename ${sp1_sample%.*})
sp2_out=${out_dir}/$(basename ${sp2_sample%.*})

# EXECUTION

# create non existent directories
mkdir -p ${out_dir}

# create popualation alle frequencies wih vcftools
vcftools --gzvcf ${ann_vcf} --keep ${sp1_sample} --freq --out ${sp1_out}
vcftools --gzvcf ${ann_vcf} --keep ${sp2_sample} --freq --out ${sp2_out}

