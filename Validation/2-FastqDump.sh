#!/bin/bash
#SBATCH -N 2
#SBATCH -c 4
#SBATCH --mem=64000
#SBATCH --ntasks-per-node=1
#SBATCH -A dutta_tumor
#SBATCH --partition=parallel
#SBATCH --time=71:00:00

SRR=$(cat /project/dutta_tumor/Ajay/Coverage/dbGAP/SRR.txt | awk NR==${SLURM_ARRAY_TASK_ID}) 
cd /project/dutta_tumor/Ajay/Coverage/dbGAP/SRA/dbGaP-24960
module load sratoolkit/2.9.1

fastq-dump --split-files $SRR -O /project/dutta_tumor/Ajay/Coverage/dbGAP/SRA/dbGaP-24960/fastq/$SRR

echo "DONE"

