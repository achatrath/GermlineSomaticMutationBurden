#!/bin/bas
#SBATCH --mem=5000
#SBATCH -A dutta_tumor
#SBATCH --partition=standard
#SBATCH --time=160:00:00

cd /project/dutta_tumor/Ajay/Coverage/dbGAP/SRA/dbGaP-24960

SRR=$(cat /project/dutta_tumor/Ajay/Coverage/dbGAP/SRA_Download.txt | cut -f2 | awk NR==${SLURM_ARRAY_TASK_ID})

module load sratoolkit/2.9.1
cd /project/dutta_tumor/Ajay/Coverage/dbGAP/SRA
prefetch $SRR

