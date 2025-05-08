#!/bin/bash
#
#SBATCH --job-name=snpAnno
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=0-01:00:00
#SBATCH --partition=basic,short
#SBATCH --output=/lisc/scratch/botany/frschmidt/mk_test/logs/snpEff_annotate%j.log
#SBATCH --error=/lisc/scratch/botany/frschmidt/mk_test/logs/snpEff_annotate%j.err
#SBATCH --mail-type=END

# ENVIRONMENT 
ml conda
conda activate diosel

SNPEFF_DIR=/lisc/scratch/botany/frschmidt/bin/snpEff
SNPEFF=${SNPEFF_DIR}/snpEff.jar

module list
java -jar $SNPEFF -version

# VARIABLES
WD=/lisc/scratch/botany/frschmidt/mk_test

# snpEff dirs and files
annot_version="vie.refAnnot"
SNPEFF_CONFIG=${SNPEFF_DIR}/snpEff.config 
db_dir=${WD}/data/snpEff_DB/

# Calciphila and sp PicN'ga vcf-file
data_dir=${WD}/data/variants
vcf=${data_dir}/pair_fla_spn/vieref_fla_spn.vcf.gz

# output directory and annotated vcf-file
out_dir=${WD}/raw_data/snpeff_annotation/pair_fla_spn
vcf_ann=${out_dir}/$(basename ${vcf%.vcf.gz}).conda.ann.vcf.gz

# EXECUTION
mkdir -p ${data_dir}
mkdir -p ${out_dir}

cd $out_dir

# run snpeff on specified vcf file
java -Xmx8g -jar ${SNPEFF} ${annot_version} -dataDir ${db_dir} ${vcf} | bgzip -c > ${vcf_ann}
