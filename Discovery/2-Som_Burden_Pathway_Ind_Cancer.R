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

overlapped_pathways <- function(result) {
	path = unique(result$perturbed_genes)
	list = list()
	for(i in 1:length(path)) {
		temp = result[result$perturbed_genes == path[i], ]
		new_pathway = paste(temp$pathway, collapse = ", ")
		temp$pathway = new_pathway
		temp = unique(temp)
		list[[i]] = temp
	}
	new_frame = as.data.frame(do.call("rbind", list), stringsAsFactors=FALSE)
	return(new_frame)
}

get_genes <- function(string) {
	tab = sort(table(strsplit2(unlist(strsplit(string, ", ")), "_", 2)))
	return(paste(paste(names(tab), " (", tab, ")", sep = ''), collapse = ", "))
}

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

pathways = readLines("/project/dutta_tumor/Ajay/PATH_GV_SIG/Data/Gene_Sets/c2.cp.reactome.v7.0.symbols.gmt")

cancers = sort(unique(PGV2$cancer))

list = list()

for(j in 1:length(cancers)) {

	p1 = rep(NA, length(pathways))
	coefficients1 = rep(NA, length(pathways))

	p2 = rep(NA, length(pathways))
	coefficients2 = rep(NA, length(pathways))

	p3 = rep(NA, length(pathways))
	coefficients3 = rep(NA, length(pathways))

	num = rep(NA, length(pathways))
	pathway = rep(NA, length(pathways))

	perturbed_genes = rep(NA, length(pathways))

	PGV_temp = PGV2[is.element(PGV2$cancer, cancers[j]), ]

	for(i in 1:length(pathways)) {
		genes = unlist(strsplit(pathways[[i]], '\t'))
		temp = PGV_temp[is.element(PGV_temp$HUGO_Symbol, genes), ]

		if(nrow(temp) < 5)
			next
		
		test = merge
		test$GV = 0
		test$GV[is.element(test$patient, temp$pt_barcode)] = 1
		test = test[complete.cases(test), ]
		if(sum(test$GV == 1) < 5)
			next

		linearMod <- try(lm(Overall_Per_MB ~ GV + age_at_initial_pathologic_diagnosis + gender +
				European + test[,9] + test[,10] + test[,11], data=test))  # build linear regression model on full data

		p1[i] = summary(linearMod)$coefficients[2,4]
		coefficients1[i] = summary(linearMod)$coefficients[2,1]

		linearMod <- try(lm(Nonsyn_Per_MB ~ GV + age_at_initial_pathologic_diagnosis + gender +
		European + test[,9] + test[,10] + test[,11], data=test))  # build linear regression model on full data


		p2[i] = summary(linearMod)$coefficients[2,4]
		coefficients2[i] = summary(linearMod)$coefficients[2,1]

		linearMod <- try(lm(Clonal_Nonsyn_Per_MB ~ GV + age_at_initial_pathologic_diagnosis + gender +
		European + test[,9] + test[,10] + test[,11], data=test))  # build linear regression model on full data


		p3[i] = summary(linearMod)$coefficients[2,4]
		coefficients3[i] = summary(linearMod)$coefficients[2,1]

		num[i] = sum(test$GV == 1)
		pathway[i] = genes[1]
		perturbed_genes[i] = paste(sort(paste(temp$bcr_patient_barcode, temp$HUGO_Symbol, sep = "_")), collapse = ", ")
	}

	result = as.data.frame(cbind(p1, coefficients1, p2, coefficients2, p3, coefficients3, num, pathway, perturbed_genes, cancers[j]), stringsAsFactors=FALSE)
	result = result[complete.cases(result), ]
	if(nrow(result) == 0)
		next
	result = overlapped_pathways(result)
	result$p1_adj = p.adjust(result$p1, method = "BH")
	result$p2_adj = p.adjust(result$p2, method = "BH")
	result$p3_adj = p.adjust(result$p3, method = "BH")

	print(j)
	list[[j]] = result
}

full = as.data.frame(do.call("rbind", list), stringsAsFactors=FALSE)

full[which(full$p3_adj < 0.05), ]

full$gene_list = unlist(lapply(full$perturbed_genes, FUN=get_genes))
full = full[order(full$p3_adj), ]

fwrite(full, "/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Results/Discovery/ind_cancer_pathway.txt", sep = '\t')

# now make plot #

data = full
colors = fread3("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/TCGA_FILES/colors.tsv", sep = '\t')
col = unlist(strsplit(colors$V1[27], ";"))


library(circlize)
library(gtools)
library(ggplot2)
library(survival)
library(survminer)

temp = data[,c(5,6,10,13)]
temp$h = -log(as.numeric(temp$p3_adj))

set.seed(1)
rows <- sample(nrow(temp))
temp = temp[rows, ]

temp = temp[order(temp$V10), ]

temp$number = 1:nrow(temp)

temp$factor = as.factor(temp$V10)

colnames(temp)[3] = "Cancer"

intercept = -log(0.05)

pdf("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Results/Discovery/Ind_cancer_manhattan.pdf", height = 7, width = 7)
p <- ggplot(temp, aes(x=temp$number, y=temp$h, color = Cancer)) + geom_point(size = 3) + scale_color_manual(values = col) + 
	geom_hline(yintercept = intercept, lwd = 1.2) + xlab("Pathway") + theme_bw(base_size=18) + ylab("-ln(FDR)") + theme(axis.ticks.x=element_blank(), axis.text.x = element_blank())
p
dev.off()

pdf("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Results/Discovery/Ind_cancer_manhatta_v2.pdf", height = 7, width = 7)
p <- ggplot(temp, aes(x=temp$number, y=temp$h, color = Cancer)) + geom_point(size = 3) + scale_color_manual(values = col) + 
	geom_hline(yintercept = intercept, lwd = 1.2) + xlab("Reactome Pathway") + theme_bw(base_size=18) + ylab("-ln(FDR)") + theme(axis.ticks.x=element_blank(), axis.text.x = element_blank())
p
dev.off()



















