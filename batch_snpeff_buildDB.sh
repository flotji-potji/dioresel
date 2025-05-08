#!/bin/bash
#
#SBATCH --job-name=snpDB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=0-00:05:00
#SBATCH --partition=basic,short
#SBATCH --output=/lisc/scratch/botany/frschmidt/mk_test/logs/snpEff_buildDB%j.log
#SBATCH --error=/lisc/scratch/botany/frschmidt/mk_test/logs/snpEff_buildDB%j.err
#SBATCH --mail-type=END

# ENVIRONMENT 
ml java
ml gffread
SNPEFF_DIR=/lisc/scratch/botany/frschmidt/bin/snpEff
SNPEFF=${SNPEFF_DIR}/snpEff.jar

module list
java -jar $SNPEFF -version

# VARIABLES
WD=/lisc/scratch/botany/frschmidt/mk_test

# snpEff dirs and files
annot_version="vie.refAnnot"
SNPEFF_CONFIG=${SNPEFF_DIR}/snpEff.config 

# vieillardii reference fasta and annotation
REF_DIR=${WD}/data/reference
REF_FASTA=${REF_DIR}/vie_reference_genome/vieillardii1167c.asm.bp.p_ctg.fa
REF_ANNOT=${REF_DIR}/vieillardii_braker_tk.gff 

# output directories and files
db_dir=${WD}/data/snpEff_DB/${annot_version}

# CDS file
cds_fa=${WD}/raw_data/transdecoder_single_best/vie_exon.fa.transdecoder.trimmed.pep.fa

# EXECUTION
mkdir -p ${db_dir}

# test if CDS fasta file already exists if not create one
if [ ! -e ${cds_fa} ]; then
	gffread -g ${REF_FASTA} -x ${cds_fa} <(grep CDS ${REF_ANNOT}) 
fi

# add new species to snpEff configuration file
grep ${annot_version} ${SNPEFF_CONFIG} || cat >> ${SNPEFF_CONFIG} <<EOF
# Diospyros vieillardii genome, version ${annot_version}
${annot_version}.genome: Diospyros vieillardii (Ebenaceae)
EOF

# add vie reference fasta sequence to database dir
ln -f -s ${REF_FASTA} ${db_dir}/sequences.fa
# add vie reference gff file to database dir
ln -f -s ${REF_ANNOT} ${db_dir}/genes.gff
# add CDS fasta to database dir
ln -f -s ${cds_fa} ${db_dir}/protein.fa

# create new vie database
cd ${db_dir}
java -jar $SNPEFF build -gff3 -v -noCheckCds -dataDir ${db_dir}/../ ${annot_version}
