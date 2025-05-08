#!/bin/bash
#
#SBATCH --job-name=blastn
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=8G
#SBATCH --time=0-00:10:00
#SBATCH --partition=basic,short
#SBATCH --output=/lisc/scratch/botany/frschmidt/mk_test/logs/blastn_genes%j.log
#SBATCH --error=/lisc/scratch/botany/frschmidt/mk_test/logs/blastn_genes%j.err
#SBATCH --mail-type=END

### ENVIRONMENT
ml ncbiblastplus
ml bedtools

### VARIABLES
wd=/lisc/scratch/botany/frschmidt/mk_test
db_dir=${wd}/data/blastdb

ref_gff=${wd}/data/reference/vieillardii_braker_tk.gff 
ref_fa=${wd}/data/reference/vie_reference_genome/vieillardii1167c.asm.bp.p_ctg.fa

raw_dir=${wd}/raw_data/reference
gene="g8348"
ident_thresh="50"
gene_fa=${raw_dir}/${gene}.fa

out_dir=${wd}/results/blastn_hits
out_hit=${out_dir}/${gene}

### EXECUTION
mkdir -p ${out_dir}

if [ ! -e ${gene_fa} ]; then
		bedtools getfasta -bed <(grep ${gene} ${ref_gff} | grep -v "Parent") \
				  		  -fi ${ref_fa} \
				  		  -fo ${gene_fa}
fi

for db in ${db_dir}/*; do
	blastn -query ${gene_fa} -db ${db}/$(basename ${db}) \
		   -out ${out_hit}.${ident_thresh}.$(basename ${db}).tsv -outfmt 6 \
		   -perc_identity ${ident_thresh}
done
