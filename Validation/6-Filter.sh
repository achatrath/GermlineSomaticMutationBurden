#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=50000
#SBATCH --ntasks-per-node=1
#SBATCH -A dutta_tumor
#SBATCH --partition=standard
#SBATCH --time=160:00:00

cd /nv/vol169/cphg_ratan/Ajay/CURRENT_OUTPUT_2020_04_03/out_dir_combined
module load bedtools

zcat /nv/vol169/cphg_ratan/Ajay/CURRENT_OUTPUT_2020_04_03/out_dir_combined/final.vcf.gz | bedtools intersect -a stdin -b /nv/vol169/cphg_ratan/Ajay/CURRENT_OUTPUT_2020_04_03/out_dir_combined/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list -sorted > output_filtered.vcf
cat output_filtered.vcf | grep PASS > output_filtered_PASS.vcf
