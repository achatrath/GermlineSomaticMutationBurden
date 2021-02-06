#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=16000
#SBATCH --ntasks-per-node=1
#SBATCH -A dutta_tumor
#SBATCH --partition=standard
#SBATCH --time=160:00:00


java -Xmx4g -Dconfig.file=/project/dutta_tumor/Ajay/Workflow/TEST_RUN/rivanna.conf -jar /project/dutta_tumor/Ajay/Workflow/cromwell-48.jar run /project/dutta_tumor/Ajay/Workflow/Workflows/b37_gvcf_to_variants.wdl --inputs /project/dutta_tumor/Ajay/Coverage/dbGAP/5-b37_gvcf_to_variants.json



