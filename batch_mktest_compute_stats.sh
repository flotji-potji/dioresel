#!/bin/bash
#
#SBATCH --job-name=mk_sum
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=0-03:00:00
#SBATCH --partition=basic,short
#SBATCH --output=/lisc/scratch/botany/frschmidt/mk_test/logs/mktest_summary%j.log
#SBATCH --error=/lisc/scratch/botany/frschmidt/mk_test/logs/mktest_summary%j.err
#SBATCH --mail-type=END

### ENVIRONMENT
ml R

module list

### VARIABLES
SCRIPTS=/lisc/project/botany/frschmidt/mk_test/scripts
WD=/lisc/scratch/botany/frschmidt/mk_test
vcf=${WD}/raw_data/snpeff_annotation/pair_fla_spn/vieref_fla_spn.ann.vcf.ann.miss-syno.vcf.gz
outgroup=${WD}/raw_data/pop_allel_freq/pair_fla_spn/vieref_fla_spn.ann.vcf.fla.frq
target=${WD}/raw_data/pop_allel_freq/pair_fla_spn/vieref_fla_spn.ann.vcf.spn.frq

### EXECUTION
R --vanilla -s -f $SCRIPTS/batch_mktest_compute_stats.R --args ${vcf} ${outgroup} ${target}
