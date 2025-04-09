#!/bin/bash
#
#SBATCH --job-name=blastmakedb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=8G
#SBATCH --time=0-00:10:00
#SBATCH --partition=basic,short
#SBATCH --output=/lisc/scratch/botany/frschmidt/mk_test/logs/makeblastdb_dioasm%j.log
#SBATCH --error=/lisc/scratch/botany/frschmidt/mk_test/logs/makeblastdb_dioasm%j.err
#SBATCH --mail-type=END

### ENVIRONMENT
ml ncbiblastplus

### VARIABLES
wd=/lisc/scratch/botany/frschmidt/mk_test
assembly_dir=/lisc/scratch/botany/ovidiu/diospyros/revio/assemblies
data_dir=${wd}/data/assemblies
out_dir=${wd}/data/blastdb

### EXECUTION
mkdir -p ${data_dir}
mkdir -p ${out_dir}

for assembly in ${assembly_dir}/*.fa; do
	asm_short=$(basename ${assembly%*.})
	asm_dir=${data_dir}/${asm_short}
	db_dir=${out_dir}/${asm_short}

	mkdir ${asm_dir}	
	mkdir ${db_dir}

	ln -f -s ${assembly} ${asm_dir}
	
	makeblastdb -in ${asm_dir}/*.fa -out ${db_dir}/${asm_short} -dbtype nucl
done
