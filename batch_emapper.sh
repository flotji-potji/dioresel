#!/bin/bash
#
#SBATCH --job-name=emapper
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=4G
#SBATCH --time=0-10:00:00
#SBATCH --partition=basic,short
#SBATCH --output=/lisc/scratch/botany/frschmidt/mk_test/logs/emapper%j.%a.log
#SBATCH --error=/lisc/scratch/botany/frschmidt/mk_test/logs/emapper%j.%a.err
#SBATCH --mail-type=END

### ENVIRONMENT
ml conda 
module list

echo "conda activate eggnog"
conda activate eggnog
emapper.py -v

echo "export EGGNOG_DATA_DIR=/lisc/scratch/mirror/eggnog-mapper/2.1.12"
export EGGNOG_DATA_DIR=/lisc/scratch/mirror/eggnog-mapper/2.1.12

### VARIABLES
wd=/lisc/scratch/botany/frschmidt/mk_test
query=${wd}/raw_data/transdecoder_single_best/chunks/chunk_${SLURM_ARRAY_TASK_ID}.pep.fa
outDir=${wd}/results/emapper/chunks
outFile=${outDir}/chunk_${SLURM_ARRAY_TASK_ID}.eggnog

### EXECUTION
echo "Started at `date`"

mkdir -p ${outDir}

emapper.py --temp_dir $TMPDIR \
           --output_dir $outDir \
           --cpu $SLURM_CPUS_PER_TASK \
           -i $query \
           -m diamond \
           --pident 60 \
           --query_cover 60 \
           --subject_cover 60 \
           --output $outFile

echo "Ended at `date`"
