#!/bin/bash
#
#SBATCH --job-name=fst_window
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --mail-type=END
#SBATCH --time=0-1:00:00
#SBATCH --output=/lisc/scratch/botany/frschmidt/fst/logs/vcftools_window%j.log
#SBATCH --error=/lisc/scratch/botany/frschmidt/fst/logs/vcftools_window%j.err

### ENVIRONMENT
ml vcftools

### VARIABLES

WD=/lisc/scratch/botany/frschmidt/fst

# directories of stored files and output directories
raw_data_dir=${WD}/raw_data
pair_name="pair_fla_umb"
pair_var_dir=${WD}/data/variants/${pair_name}
out_dir=${WD}/results/vcftools/${pair_name}

# files used for execution
vcf=${pair_var_dir}/vieref_fla_umb.vcf.gz
sp1_sample=${pair_var_dir}/flavocarpa.samples
sp2_sample=${pair_var_dir}/umbrosa.samples

out_file=${out_dir}/$(basename ${vcf%.*}).sliding.window

# variables to define windows
window_size=10000
window_step=5000

### EXECUTION
mkdir -p ${out_dir}

vcftools --gzvcf ${vcf} \
		--weir-fst-pop ${sp1_sample} \
		--weir-fst-pop ${sp2_sample} \
		--fst-window-size ${window_size} \
		--fst-window-step ${window_step} \
		--out ${out_file}
