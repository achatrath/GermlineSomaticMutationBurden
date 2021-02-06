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

MAF_global_driver <- function(MAF_driver, CDR2, drivers) {
	global = rep(NA, nrow(CDR2))
	for(i in 1:nrow(CDR2)) {
		print(i)
		type = CDR2$type[i]
		if(type == "COAD" | type == "READ")
			type = "COADREAD"
		type = c(type, "PANCAN")

		drivers_temp = drivers[is.element(drivers$Cancer, type), ]
		temp = MAF_driver[is.element(MAF_driver$V1, drivers_temp$Gene), ]
		temp = temp[temp$patient == CDR2$pt_barcode[i], ]
		if(nrow(temp) == 0)
			global[i] = 0
		else if(nrow(temp) > 0) {
			global[i] = nrow(temp)
		}
	}
	return(global)
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


MAF_driver_specific <- function(MAF_driver, gene, merge1) {
	temp = MAF_driver[MAF_driver$V1 == gene, ]
	merge1$SOM = as.numeric(is.element(merge1$pt_barcode, temp$patient))
	return(merge1)
}

getPathways <- function(pathway) {
  return(unlist(strsplit(pathway, "\t"))[1])
}

getGenes <- function(pathway) {
  return(unlist(strsplit(pathway, "\t")))
}

getBarCode <- function(code) {
  return(paste(unlist(strsplit(code, "-"))[1:3], collapse = "-"))
}

checkPathways <- function(pathways, gene) {
	found = c()

	for(i in 1:length(pathways)) {
		temp = unlist(strsplit(pathways[[i]], '\t'))
		if(is.element(gene, temp))
			found = c(found, i)
	}

	return(found)
}

getCounts <- function(merge1, temp2) {
	merge1$PATH_SOM = 0
	for(i in 1:nrow(merge1)) {
		temporary = temp2[temp2$pt_barcode == merge1$bcr_patient_barcode[i], ]
		if(nrow(temporary) == 0)
			next
		else {
			merge1$PATH_SOM[i] = nrow(temporary)
		}
	}

	return(merge1)
}


setwd("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/TCGA_FILES")

CDR = fread3('/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/TCGA_FILES/TCGA-CDR.csv', sep = '\t')
CDR2 = CDR[, c(2:5)]
CDR2$patient = strsplit2(CDR2[,1], "-", 3)

SOM2 = fread3("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Somatic_Burden_Norm_Cov.tsv", sep = '\t')
TCGAA = fread3("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/TCGA_FILES/STRUCTURE.txt", sep = '\t')
PGV = fread3("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/TCGA_FILES/PCA_pathVar_integrated_filtered_adjusted.tsv", sep = '\t')
PGV2 = PGV[,c(1:3)]
PGV2$pt_barcode = strsplit2(PGV2[,1], "-", 3)
PGV2$path_variant = 1

SOM2$patient = strsplit2(SOM2$Barcode, "-", 3)
TCGAA$patient = strsplit2(TCGAA[,1], "-", 3)

merge = merge(CDR2, TCGAA, by="patient")
merge = merge(merge, SOM2, by = "patient")

MAF = fread3("/project/dutta_tumor/Ajay/Coverage/data/mc3.v0.2.8.PUBLIC.maf", sep = '\t')

MAF_damage = unlist(lapply(MAF$PolyPhen, FUN=grepl, pattern = "damaging"))
MAF_damaged = MAF[MAF_damage, ]
MAF_damaged$patient = strsplit2(MAF_damaged$Tumor_Sample_Barcode, "-", 3)

PGV2 = PGV[,c(1:3)]
PGV2$pt_barcode = strsplit2(PGV2[,1], "-", 3)
PGV2$path_variant = 1

CDR2 = CDR[,c(2,5)]
CDR2$pt_barcode = strsplit2(CDR2[,1], "-", 3)

tab = sort(table(MAF_damaged$Hugo_Symbol), decreasing = TRUE)
genes = names(tab)

temp = PGV2[PGV2$HUGO_Symbol == "PMS2", ]
merge1 = merge
merge1$PATH_GV = 0
merge1$PATH_GV[is.element(merge1$patient, temp$pt_barcode)] = 1

p = rep(NA, length(genes))
coefficients = rep(NA, length(genes))
mean1= rep(NA, length(genes))
mean2 = rep(NA, length(genes))

library(doParallel)
library(doSNOW)
cl <- makeCluster(9, outfile = "")
registerDoParallel(cl)

list <- foreach(i=1:length(genes)) %dopar% {
	if(i%%100 == 0)
		print(i)
	temp2 = MAF_damaged[MAF_damaged$Hugo_Symbol == genes[i], ]
	temp2 = temp2[strsplit2(temp2[,16], "-", 4) == "01A" | strsplit2(temp2[,16], "-", 4) == "03A", ]

	if(nrow(temp2) == 0)
		return(NA)
	merge1$PATH_SOM = 0
	merge1$PATH_SOM[is.element(merge1$patient, temp2$patient)] = 1
	merge1$PATH_SOM = factor(merge1$PATH_SOM)

	model <- glm(PATH_SOM ~ PATH_GV + type + Clonal_Nonsyn_Per_MB + 
		age_at_initial_pathologic_diagnosis + gender + merge1[,8] + merge1[,9] +
		merge1[,10] + merge1[,11], data = merge1, family = "binomial")  
	p = summary(model)$coefficients[2,4]
	coefficients = summary(model)$coefficients[2,1]
	mean1= sum(as.numeric(merge1$PATH_SOM[merge1$PATH_GV == 1])-1)/length(merge1$PATH_SOM[merge1$PATH_GV == 1])
	mean2= sum(as.numeric(merge1$PATH_SOM[merge1$PATH_GV == 0])-1)/length(merge1$PATH_SOM[merge1$PATH_GV == 0])
	num = sum(as.numeric(merge1$PATH_SOM))
	res = as.data.frame(cbind(genes[i], coefficients, p, mean1, mean2, num), stringsAsFactors=FALSE)
	return(res)
}

full = as.data.frame(do.call("rbind", list), stringsAsFactors=FALSE)

fwrite(full, "/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Results/Molecular/Control/PMS2_other_genes.txt", sep = '\t')


