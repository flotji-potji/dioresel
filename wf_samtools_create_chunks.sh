#!/bin/bash

ml samtools

list=($(cut -f 1 $2))
num_trans=$(wc -l < $2)
size_chunk=1000
num_chunks=$(($num_trans/$size_chunk))

for i in $(seq 0 "$num_chunks"); do
	si=$(($i * "$size_chunk"))
	samtools faidx $1 ${list[@]:$si:$size_chunk} > chunk_$(($i+1)).pep.fa
done
