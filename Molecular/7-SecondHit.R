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

genes = c("APC", "FANCL", "SLC25A13", "ERCC3", "MSH6", "PMS2", "TP53")
p = rep(NA, length(genes))
coefficients = rep(NA, length(genes))
mean1= rep(NA, length(genes))
mean2 = rep(NA, length(genes))



for(i in 1:length(genes)) {
	print(i)
	temp = PGV2[PGV2$HUGO_Symbol == genes[i], ]
	merge1 = merge
	merge1$PATH_GV = 0
	merge1$PATH_GV[is.element(merge1$patient, temp$pt_barcode)] = 1
	temp2 = MAF_damaged[MAF_damaged$Hugo_Symbol == genes[i], ]
	temp2 = temp2[strsplit2(temp2[,16], "-", 4) == "01A" | strsplit2(temp2[,16], "-", 4) == "03A", ]

	if(nrow(temp2) == 0)
		next
	merge1$PATH_SOM = 0
	merge1$PATH_SOM[is.element(merge1$patient, temp2$patient)] = 1
	merge1$PATH_SOM = factor(merge1$PATH_SOM)

	model <- glm(PATH_SOM ~ PATH_GV + type + Clonal_Nonsyn_Per_MB + 
		age_at_initial_pathologic_diagnosis + gender + merge1[,8] + merge1[,9] +
		merge1[,10] + merge1[,11], data = merge1, family = "binomial")  
	p[i] = summary(model)$coefficients[2,4]
	coefficients[i] = summary(model)$coefficients[2,1]
	mean1[i] = sum(as.numeric(merge1$PATH_SOM[merge1$PATH_GV == 1])-1)/length(merge1$PATH_SOM[merge1$PATH_GV == 1])
	mean2[i] = sum(as.numeric(merge1$PATH_SOM[merge1$PATH_GV == 0])-1)/length(merge1$PATH_SOM[merge1$PATH_GV == 0])
}

res = as.data.frame(cbind(genes, coefficients, p, mean1, mean2), stringsAsFactors=FALSE)
res$p_adj = p.adjust(res$p, method = "BH")

# individual genes and pathway that the gene is contained in #

genes = c("APC", "FANCL", "SLC25A13", "ERCC3", "MSH6", "PMS2", "TP53")

pathways = readLines("/project/dutta_tumor/Ajay/PATH_GV_SIG/Data/Gene_Sets/c2.cp.reactome.v7.0.symbols.gmt")
pathway_annot = unlist(lapply(pathways, FUN=getPathways))
list = list()

for(i in 1:length(genes)) {
	print(i)
	path_search = checkPathways(pathways, genes[i])
	p = rep(NA, length(path_search))
	coefficients = rep(NA, length(path_search))
	mean_GV = rep(NA, length(path_search))
	mean_no_GV = rep(NA, length(path_search))


	for(j in 1:length(path_search)) {
		print(j)
		gene_set = getGenes(pathways[[path_search[j]]])
		merge1 = merge
		temp = PGV2[PGV2$HUGO_Symbol == genes[i], ]
		merge1$PATH_GV = 0
		merge1$PATH_GV[is.element(merge1$patient, temp$pt_barcode)] = 1

		temp2 = MAF_damaged[is.element(MAF_damaged$Hugo_Symbol, gene_set), ]
		temp2$pt_barcode = unlist(lapply(temp2[,16], FUN=getBarCode))
		temp2 = temp2[strsplit2(temp2[,16], "-", 4) == "01A" | strsplit2(temp2[,16], "-", 4) == "03A", ]
		tab = table(temp2$patient)
		tab = as.data.frame(cbind(names(tab), as.numeric(tab)), stringsAsFactors=FALSE)
		colnames(tab) = c("patient", "PATH_SOM")

		merge2 = merge(merge1, tab, by = "patient", all = TRUE)
		merge2$PATH_SOM[is.na(merge2$PATH_SOM)] = 0
		merge2$PATH_SOM = as.numeric(merge2$PATH_SOM)

		model <- lm(PATH_SOM ~ PATH_GV + type + Clonal_Nonsyn_Per_MB + 
			age_at_initial_pathologic_diagnosis + gender + merge2[,8] + merge2[,9] +
			merge2[,10] + merge2[,11] + Coverage, data = merge2)

		model2 <- lm(PATH_SOM ~ PATH_GV, data = merge2)

		p[j] = summary(model)$coefficients[2,4]
		coefficients[j] = summary(model)$coefficients[2,1]
		mean_GV[j] = mean(merge2$PATH_SOM[merge1$PATH_GV == 1] / merge2$Nonsyn_Per_MB[merge1$PATH_GV == 1], na.rm = TRUE)
		mean_no_GV[j] = mean(merge2$PATH_SOM[merge1$PATH_GV == 0] / merge2$Nonsyn_Per_MB[merge1$PATH_GV == 0], na.rm = TRUE)
	}

	res = as.data.frame(cbind(genes[i], pathway_annot[path_search], p, coefficients, mean_GV, mean_no_GV), stringsAsFactors=FALSE)
	res$padj = p.adjust(res$p, method = "BH")
	list[[i]] = res
}

results = as.data.frame(do.call("rbind", list), stringsAsFactors=FALSE)
results$padj2 = p.adjust(results$p, method = "BH")
results = results[order(results$padj2), ]
results2 = results[results$padj2 < 0.05 & results$coefficients > 0, ]

table(results2$mean_GV > results2$mean_no_GV)
results2 = results2[results2$mean_GV > results2$mean_no_GV, ]

fwrite(results2, "/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Results/Molecular/Somatic_Mutations/results_gene_pathway.txt", sep = '\t')

# Pathways in individual cancers #

prev_result = fread3("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Results/Discovery/ind_cancer_pathway.txt", sep = '\t')
prev_result = prev_result[prev_result$p3_adj < 0.05, ]
reformatted_prev = reformat_prev_result(prev_result)

pathways = readLines("/project/dutta_tumor/Ajay/PATH_GV_SIG/Data/Gene_Sets/c2.cp.reactome.v7.0.symbols.gmt")
pathway_annot = unlist(lapply(pathways, FUN=getPathways))
list = list()

for(i in 1:nrow(reformatted_prev)) {
	print(i)
	path_search = reformatted_prev[i, 1]
	p = rep(NA, length(path_search))
	coefficients = rep(NA, length(path_search))
	mean_GV = rep(NA, length(path_search))
	mean_no_GV = rep(NA, length(path_search))

	for(j in 1:length(path_search)) {
		print(j)
		index = which(pathway_annot == path_search[j])
		gene_set = getGenes(pathways[[index]])
		merge1 = merge
		temp = PGV2[is.element(PGV2$HUGO_Symbol, gene_set), ]
		merge1$PATH_GV = 0
		merge1$PATH_GV[is.element(merge1$patient, temp$pt_barcode)] = 1

		temp2 = MAF_damaged[is.element(MAF_damaged$Hugo_Symbol, gene_set), ]
		temp2$pt_barcode = unlist(lapply(temp2[,16], FUN=getBarCode))
		temp2 = temp2[strsplit2(temp2[,16], "-", 4) == "01A" | strsplit2(temp2[,16], "-", 4) == "03A", ]
		tab = table(temp2$patient)
		tab = as.data.frame(cbind(names(tab), as.numeric(tab)), stringsAsFactors=FALSE)
		colnames(tab) = c("patient", "PATH_SOM")

		merge2 = merge(merge1, tab, by = "patient", all = TRUE)
		merge2$PATH_SOM[is.na(merge2$PATH_SOM)] = 0
		merge2$PATH_SOM = as.numeric(merge2$PATH_SOM)

		merge2 = merge2[merge2$type == reformatted_prev[i,2], ]

		if(reformatted_prev[i,2] == "UCEC")
			model <- lm(PATH_SOM ~ PATH_GV + Clonal_Nonsyn_Per_MB + 
				age_at_initial_pathologic_diagnosis + merge2[,8] + merge2[,9] +
				merge2[,10] + merge2[,11] + Coverage, data = merge2)
		else
			model <- lm(PATH_SOM ~ PATH_GV + Clonal_Nonsyn_Per_MB + 
				age_at_initial_pathologic_diagnosis + merge2[,8] + merge2[,9] +
				merge2[,10] + merge2[,11] + Coverage + gender, data = merge2)


		model2 <- lm(PATH_SOM ~ PATH_GV, data = merge2)

		p[j] = summary(model)$coefficients[2,4]
		coefficients[j] = summary(model)$coefficients[2,1]
		mean_GV[j] = mean(merge2$PATH_SOM[merge1$PATH_GV == 1] / merge2$Nonsyn_Per_MB[merge1$PATH_GV == 1], na.rm = TRUE)
		mean_no_GV[j] = mean(merge2$PATH_SOM[merge1$PATH_GV == 0] / merge2$Nonsyn_Per_MB[merge1$PATH_GV == 0], na.rm = TRUE)
	}

	res = as.data.frame(cbind(path_search, reformatted_prev[i,2], p, coefficients, mean_GV, mean_no_GV), stringsAsFactors=FALSE)
	list[[i]] = res
}

results = as.data.frame(do.call("rbind", list), stringsAsFactors=FALSE)
results$padj2 = p.adjust(results$p, method = "BH")
results = results[order(results$padj2), ]
results2 = results[results$padj2 < 0.05 & results$coefficients > 0, ]

table(results2$mean_GV > results2$mean_no_GV)
results2 = results2[results2$mean_GV > results2$mean_no_GV, ]

fwrite(results2, "/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Results/Molecular/Somatic_Mutations/results_pathway_ind_cancer.txt", sep = '\t')

# pancancer pathway # 

prev_result = fread3("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Results/Discovery/pancancer_pathway.txt", sep = '\t')
prev_result = prev_result[prev_result$p3_adj < 0.05, ]
reformatted_prev = reformat_prev_result(prev_result)

pathways = readLines("/project/dutta_tumor/Ajay/PATH_GV_SIG/Data/Gene_Sets/c2.cp.reactome.v7.0.symbols.gmt")
pathway_annot = unlist(lapply(pathways, FUN=getPathways))
list = list()

for(i in 1:nrow(reformatted_prev)) {
	print(i)
	path_search = reformatted_prev[i, 1]
	p = rep(NA, length(path_search))
	coefficients = rep(NA, length(path_search))
	mean_GV = rep(NA, length(path_search))
	mean_no_GV = rep(NA, length(path_search))

	for(j in 1:length(path_search)) {
		print(j)
		index = which(pathway_annot == path_search[j])
		gene_set = getGenes(pathways[[index]])
		merge1 = merge
		temp = PGV2[is.element(PGV2$HUGO_Symbol, gene_set), ]
		merge1$PATH_GV = 0
		merge1$PATH_GV[is.element(merge1$patient, temp$pt_barcode)] = 1

		temp2 = MAF_damaged[is.element(MAF_damaged$Hugo_Symbol, gene_set), ]
		temp2$pt_barcode = unlist(lapply(temp2[,16], FUN=getBarCode))
		temp2 = temp2[strsplit2(temp2[,16], "-", 4) == "01A" | strsplit2(temp2[,16], "-", 4) == "03A", ]
		tab = table(temp2$patient)
		tab = as.data.frame(cbind(names(tab), as.numeric(tab)), stringsAsFactors=FALSE)
		colnames(tab) = c("patient", "PATH_SOM")

		merge2 = merge(merge1, tab, by = "patient", all = TRUE)
		merge2$PATH_SOM[is.na(merge2$PATH_SOM)] = 0
		merge2$PATH_SOM = as.numeric(merge2$PATH_SOM)

		model <- lm(PATH_SOM ~ PATH_GV + Clonal_Nonsyn_Per_MB + 
			age_at_initial_pathologic_diagnosis + merge2[,8] + merge2[,9] +
			merge2[,10] + merge2[,11] + Coverage + gender, data = merge2)


		model2 <- lm(PATH_SOM ~ PATH_GV, data = merge2)

		p[j] = summary(model)$coefficients[2,4]
		coefficients[j] = summary(model)$coefficients[2,1]
		mean_GV[j] = mean(merge2$PATH_SOM[merge1$PATH_GV == 1] / merge2$Nonsyn_Per_MB[merge1$PATH_GV == 1], na.rm = TRUE)
		mean_no_GV[j] = mean(merge2$PATH_SOM[merge1$PATH_GV == 0] / merge2$Nonsyn_Per_MB[merge1$PATH_GV == 0], na.rm = TRUE)
	}

	res = as.data.frame(cbind(path_search, p, coefficients, mean_GV, mean_no_GV), stringsAsFactors=FALSE)
	list[[i]] = res
}

results = as.data.frame(do.call("rbind", list), stringsAsFactors=FALSE)
results$padj2 = p.adjust(results$p, method = "BH")
results = results[order(results$padj2), ]
results2 = results[results$padj2 < 0.05 & results$coefficients > 0, ]

table(results2$mean_GV > results2$mean_no_GV)
results2 = results2[results2$mean_GV > results2$mean_no_GV, ]

fwrite(results2, "/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Results/Molecular/Somatic_Mutations/results_pathway_ind_cancer.txt", sep = '\t')






