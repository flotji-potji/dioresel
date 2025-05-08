#!/bin/bash
#
#SBATCH --job-name=agatME
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=8G
#SBATCH --time=0-01:00:00
#SBATCH --partition=basic,short
#SBATCH --output=/lisc/scratch/botany/frschmidt/mk_test/logs/agat_extract_exons%j.log
#SBATCH --error=/lisc/scratch/botany/frschmidt/mk_test/logs/agat_extract_exons%j.err
#SBATCH --mail-type=END

### ENVIRONMENT
ml conda
conda activate agat

module list
echo "agat version $(agat -v)"

### VARIABLES
WD=/lisc/scratch/botany/frschmidt/mk_test

# vieillardii reference fasta and annotation
REF_DIR=${WD}/data/reference
REF_FASTA=${REF_DIR}/vie_reference_genome/vieillardii1167c.asm.bp.p_ctg.fa
REF_ANNOT=${REF_DIR}/vieillardii_braker_tk.gff 

# output directories and files
tmp_fasta=${TMPDIR}/$(basename ${REF_FASTA})

raw_out_dir=${WD}/raw_data/reference
exon_fa=${raw_out_dir}/vie_exon.fa # fasta-file of vie genes expanded 1000bp


### EXECUTION
mkdir -p ${raw_out_dir}

# create a multiline fasta of reference genome
# also serves as plase for agat to store index file into
fold ${REF_FASTA} > ${tmp_fasta}

cd ${raw_out_dir}

agat_sp_extract_sequences.pl --gff ${REF_ANNOT} \
							--fasta ${tmp_fasta} \
							-t exon --merge -o ${exon_fa}

