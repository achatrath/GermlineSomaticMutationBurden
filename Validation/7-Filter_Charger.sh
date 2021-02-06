

cat charger_summary.tsv | cut -f2,3,4,20 | grep Pathogenic > pathogenic_charger.tsv
cat charger_summary.tsv | cut -f2,3,4,20,1 | grep Pathogenic > pathogenic_gene_charger.tsv
cd /project/dutta_tumor/Ajay/Coverage/dbGAP/Workflow/out_dir_combined
module load gcc
module load bedtools
zcat final.vcf.gz | bedtools intersect -a stdin -b /project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/pathogenic_charger.tsv > /project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/pathogenic.vcf