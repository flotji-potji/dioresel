#!/bin/bash
#
#SBATCH --job-name=interpro
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=0-04:00:00
#SBATCH --partition=basic,short
#SBATCH --output=/lisc/scratch/botany/frschmidt/mk_test/logs/interpro%j.%a.log
#SBATCH --error=/lisc/scratch/botany/frschmidt/mk_test/logs/interpro%j.%a.err
#SBATCH --mail-type=END
#SBATCH --license=interpro

### ENVIRONMENT
ml interproscan
module list

### VARIABLES
wd=/lisc/scratch/botany/frschmidt/mk_test

query=${wd}/raw_data/transdecoder_single_best/chunks/chunk_${SLURM_ARRAY_TASK_ID}.pep.fa
outDir=${wd}/results/interpro/chunks
outFile=${outDir}/chunk_${SLURM_ARRAY_TASK_ID}.tsv

### EXECUTION
echo "Started at `date`"

mkdir -p ${outDir}
cd ${outDir}

cmd=$(cat <<EOF
interproscan.sh -i ${query} -b ${outFile} -cpu $SLURM_CPUS_PER_TASK -f TSV -goterms -pa -iprlookup -T ${TMPDIR} -t p
EOF
)
echo $cmd
/usr/bin/time -format="%E\tTIME\t%M\tMEM\t%P\tCPU\t$cmd" $cmd

echo "Ended at `date`"
