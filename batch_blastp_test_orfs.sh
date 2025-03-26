#!/bin/bash
#
#SBATCH --job-name=blastp
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=4G
#SBATCH --time=0-01:00:00
#SBATCH --partition=basic,short
#SBATCH --output=/lisc/scratch/botany/frschmidt/mk_test/logs/blastp_test_orfs%j.%a.log
#SBATCH --error=/lisc/scratch/botany/frschmidt/mk_test/logs/blastp_test_orfs%j.%a.err
#SBATCH --mail-type=END
#SBATCH --license=scratch-highio

### ENVIRONMENT
ml ncbiblastplus
module list

### VARIABLES
WD=/lisc/scratch/botany/frschmidt/mk_test
db=/lisc/scratch/mirror/ncbi/2025-02-25/nr
trans=${WD}/results/complete_orf_transdecoder/vie_exon.fa.transdecoder_dir/chunks/chunk_${SLURM_ARRAY_TASK_ID}.pep.fa
outDir=${WD}/raw_data/blast_hits/chunks
outFile=${outDir}/chunk_${SLURM_ARRAY_TASK_ID}.blastp.outfmt6

### EXECUTION
echo "Started at `date`"

mkdir -p $outDir
echo "mkdir -p $outDir"
cd $outDir
echo "cd $outDir"

blastp -query $trans -db $db -out $outFile -num_threads $SLURM_CPUS_PER_TASK -outfmt 6 -max_target_seqs 1 -evalue 1e-5

echo "Ended at `date`"
