fread2 <- function(name, sep) {
  library(data.table)
  data = fread(name, sep = sep)
  data = as.data.frame(data, stringsAsFactors=FALSE)
  rownames(data) = data[,1]
  data = data[,-1]
  return(data)
}

fread3 <- function(name, sep) {
  library(data.table)
    data = fread(name, sep = sep)
    data = as.data.frame(data, stringsAsFactors=FALSE)
    return(data)
}

strsplit2 <- function(vector, character, position) {
        values = rep('a', length(vector))
        for(i in 1:length(vector))
                values[i] = unlist(strsplit(vector[i], character))[position]
        return(values)
}

fastsplit <- function(string) {
  return(paste(unlist(strsplit(string, "-"))[1:3], collapse = "-"))
}

get_pathways_list <- function(pathways) {
  vec <- rep(NA, length(pathways))
  for(i in 1:length(pathways)) {
    vec[i] = unlist(strsplit(pathways[[i]], "\t"))[1]
  }
  return(vec)
}

reformat_prev_result <- function(prev_result) {
  new = prev_result[,c(8,10)]
  colnames(new) = c("pathway", "cancer")
  list = list()

  for(i in 1:nrow(new)) {
    temp = unlist(strsplit(new$pathway[i], ", "))
    new_temp = as.data.frame(cbind(temp, new$cancer[i]), stringsAsFactors=FALSE)
    colnames(new_temp) = c("pathway", "cancer")
    list[[i]] = new_temp
  }

  frame = as.data.frame(do.call("rbind", list), stringsAsFactors=FALSE)
  return(frame)
}

library(signeR)
library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg19)
library(maftools)
library(BSgenome)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
library(MutationalPatterns)

MAF = fread3("/project/dutta_tumor/Ajay/Coverage/data/mc3.v0.2.8.PUBLIC.maf", sep = '\t')
MAF$patient = unlist(mclapply(MAF[,16], FUN=fastsplit, mc.cores = 8))

setwd("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/TCGA_FILES")

CDR = fread3('/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/TCGA_FILES/TCGA-CDR.csv', sep = '\t')
CDR2 = CDR[, c(2:5)]
CDR2$patient = strsplit2(CDR2[,1], "-", 3)

TCGAA = fread3("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/TCGA_FILES/STRUCTURE.txt", sep = '\t')
PGV = fread3("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/TCGA_FILES/PCA_pathVar_integrated_filtered_adjusted.tsv", sep = '\t')
PGV2 = PGV[,c(1:3)]
PGV2$pt_barcode = strsplit2(PGV2[,1], "-", 3)
PGV2$path_variant = 1

TCGAA$patient = strsplit2(TCGAA[,1], "-", 3)

merge = merge(CDR2, TCGAA, by="patient")


PGV = fread3("PCA_pathVar_integrated_filtered_adjusted.tsv", sep = '\t')
CDR = fread3('TCGA-CDR.csv', sep = '\t')

PGV2 = PGV[,c(1:3)]
PGV2$pt_barcode = strsplit2(PGV2[,1], "-", 3)

gene_names = c("APC", "FANCL", "SLC25A13", "ERCC3", "MSH6", "TP53", "PMS2", "MSH2", "BRIP1")

MAF_GV = MAF
setwd("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/TCGA_FILES")
fwrite(MAF_GV, "all_pt.maf", sep = '\t')
library(maftools)
GV = read.maf(maf = "all_pt.maf")
GV.tnm = trinucleotideMatrix(maf = GV, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
GV.tnm1 = t(GV.tnm[[1]])
mut_mat = GV.tnm1
mut_mat = mut_mat + 0.001
mut_mat = as.matrix(mut_mat)

cancer_signatures = fread3("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/TCGA_FILES/sigProfiler_exome_SBS_signatures.csv", sep = ',')
remove = c(50, 52:67)
cancer_signatures = cancer_signatures[, -remove]
A = substring(cancer_signatures$SubType, 1, 1)
B = substring(cancer_signatures$SubType, 3, 3)
C = paste(A, "[", cancer_signatures$Type, "]", B, sep = '')
rownames(cancer_signatures) = C
cancer_signatures = cancer_signatures[, -c(1:2)]
cancer_signatures = as.matrix(cancer_signatures)

# sp_url <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/","signatures_probabilities.txt", sep = "")
# cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
# new_order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
# cancer_signatures = cancer_signatures[as.vector(new_order),]
# row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
# cancer_signatures = as.matrix(cancer_signatures[,4:33])

fit_res <- fit_to_signatures(mut_mat, cancer_signatures)

result = fit_res$contribution
result2 = t(result)
sums = rowSums(result2)
result3 = result2 / sums
result3 = as.data.frame(result3, stringsAsFactors=FALSE)
result3$bcr_patient_barcode = unlist(lapply(rownames(result3), FUN=fastsplit))

# Pancancer Pathway #
prev_result = fread3("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Results/Discovery/pancancer_pathway.txt", sep = '\t')
prev_result = prev_result[prev_result$p3_adj < 0.05, ]
reformatted_prev = reformat_prev_result(prev_result)

pathways = readLines("/project/dutta_tumor/Ajay/PATH_GV_SIG/Data/Gene_Sets/c2.cp.reactome.v7.0.symbols.gmt")
pathway_list = get_pathways_list(pathways)

keep <- c("REACTOME_MISMATCH_REPAIR", 
  "REACTOME_DISEASES_OF_MISMATCH_REPAIR_MMR")

j = 1
reformatted_prev = reformatted_prev[is.element(reformatted_prev$pathway, keep), ]
index = which(pathway_list == reformatted_prev$pathway[j])
genes = unlist(strsplit(pathways[[index]], "\t"))

PGV_temp = PGV2[is.element(PGV2$HUGO_Symbol, genes), ]

returned_results = list()

for(i in 1:48) {
  temp = result3[, c(49, i)]
  temp$pt_barcode = strsplit2(temp$bcr_patient_barcode, "-", 3)
  temp$GV = 0
  temp$GV[is.element(temp$pt_barcode, PGV_temp$pt_barcode)] = 1

  merged = merge(merge, temp, by = "bcr_patient_barcode")

  model = lm(merged[,15] ~ merged$GV + merged$age_at_initial_pathologic_diagnosis + 
      merged$gender + merged[,8] + merged[,9] + merged[,10] + merged[,11])

  fold = mean(temp[temp$GV==1,2]) /  mean(temp[temp$GV==0,2])

 returned_results[[i]] = c(summary(model)$coefficients[2, 1], summary(model)$coefficients[2, 4], fold)
}

frame = as.data.frame(do.call("rbind", returned_results), stringsAsFactors=FALSE)
colnames(frame) = c("GV_estimate","GV_p", "fold")
rownames(frame) = colnames(result3)[-ncol(result3)]

frame[which(frame$GV_p < 0.05), ]
print(reformatted_prev$pathway[j])



