#!/bin/bash

# ENVIRONMENT 
ml vcftools
ml bcftools

module list

# VARIABLES

WD=/lisc/scratch/botany/frschmidt/mk_test

# directories of stored files and output directories
raw_data_dir=${WD}/raw_data
pair_name="pair_fla_spn"
pair_var_dir=${WD}/data/variants/${pair_name}
pair_ann_dir=${raw_data_dir}/snpeff_annotation/${pair_name}
out_dir=${raw_data_dir}/pop_allel_freq/${pair_name}

# files used for execution
ann_vcf=${pair_ann_dir}/vieref_fla_spn.ann.vcf.gz
ann_ms_vcf=${pair_ann_dir}/$(basename ${ann_vcf%.*}).ann.miss-syno.vcf.gz
sp1_sample=${pair_var_dir}/fla.samples
sp2_sample=${pair_var_dir}/spn.samples
sp1_out=${out_dir}/$(basename ${ann_vcf%.*}).$(basename ${sp1_sample%.*})
sp2_out=${out_dir}/$(basename ${ann_vcf%.*}).$(basename ${sp2_sample%.*})

# EXECUTION

# create non existent directories
mkdir -p ${out_dir}

# subset vcf file to only include missense and synonymous SNPs
bcftools view -i 'ANN ~ "missense" || ANN ~ "synonymous"' ${ann_vcf} -Oz -o ${ann_ms_vcf}

# create popualation alle frequencies wih vcftools
vcftools --gzvcf ${ann_ms_vcf} --keep ${sp1_sample} --freq --out ${sp1_out}
vcftools --gzvcf ${ann_ms_vcf} --keep ${sp2_sample} --freq --out ${sp2_out}

