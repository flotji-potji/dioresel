# download genes as fasta files
ml edirect
for i in {3..52}; do 
	efetch -db nucleotide -format fasta -id <(printf '%s%02d' "KT5872" $i) >> dio_genes.fasta
done
# blast genes and transform results into bed-file
blastn -query dio_genes.fa -db data/blastdb/vieillardii1167c.asm.bp.p_ctg.fa/vieillardii1167c.asm.bp.p_ctg.fa -outfmt 6 -evalue 1e-150 | cut -f2,9,10,1 > matched.bed
awk -v OFS='\t' '{if($4 < $3) print $2, $4, $3, $1; else print $2, $3, $4, $1}' matched.bed > matched_right.bed
# annotate bed file to genes
bedtools intersect -b <(awk '$3 == "gene"' ../dioresel/data/reference/vieillardii.gff) -a matched_right.bed -wo | cut -f1,2,3,4,13 | awk -v OFS='\t' '{split($5, a, "="); print $1, $2, $3, $4, a[2]}' > parallel_genes.bed

