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

var_splitter <- function(string) {
	split1=strsplit2(string, ":", 1)
	split2 = unlist(strsplit(split1, "/"))
	split3 = unlist(strsplit(split1, "[|]"))
	if((split2[1] == "." & split2[2] == ".") | (split3[1] == "." & split3[2] == "."))
		return(0)
	val1 = sum(as.numeric(unlist(strsplit(split1, "/"))))
	val2 = sum(as.numeric(unlist(strsplit(split1, "[|]"))))
	return(max(val1, val2, na.rm = TRUE))
}

not_zero <- function(vals) {
	return(sum(vals != 0))
}

not_zero_v2 <- function(vals) {
	return(max(as.numeric(vals != 0)))
}

get_pathways <- function(prev_result) {
	list = list()
	for(i in 1:nrow(prev_result)) {
		list[[i]] = unlist(strsplit(prev_result$pathway[i], ", "))
	}
	full = unlist(list)
	return(full)
}

get_genes <- function(pathways, full) {
	list = list()

	for(i in 1:length(pathways)) {
		split1 = unlist(strsplit(pathways[i], "\t"))
		if(is.element(split1[1], full))
			list[[i]] = split1
	}

	gene_list = unlist(list)
	return(gene_list)
}

VAR <- fread3("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/pathogenic.vcf", sep = '\t')
CHARGER <- fread3("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/pathogenic_gene_charger.tsv", sep = '\t')
header = fread3("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/header.txt", sep = '\t')

colnames(VAR) = colnames(header)

VAR$KEY = paste(VAR[,1], VAR[,2], sep = ".")
CHARGER$KEY = paste(CHARGER[,2], CHARGER[,3], sep = '.')

merged = merge(VAR, CHARGER, by = "KEY")

MUT = merged[,11:154]

MUT2=apply(MUT, MARGIN=c(1,2), FUN=var_splitter)

res = apply(MUT2, FUN=not_zero, MARGIN=1)
res2 = as.data.frame(cbind(res, merged$V1), stringsAsFactors=FALSE)
res2$res = as.numeric(res2$res)

res3= aggregate(res2$res ~ res2$V2, data=res2, sum)
colnames(res3) = c("Gene", "Number")
res3 = res3[order(res3$Number, decreasing = TRUE), ]

clinical = fread3("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Clinical_Melanoma_VanAllen.txt", sep = '\t')
meta = fread3("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/SRA_Download.txt", sep = '\t')

# PFS pooling pathways together #

pathways = readLines("/project/dutta_tumor/Ajay/PATH_GV_SIG/Data/Gene_Sets/c2.cp.reactome.v7.0.symbols.gmt")

prev_result = fread3("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Results/Discovery/ind_cancer_pathway.txt", sep = '\t')
prev_result = prev_result[prev_result$p3_adj < 0.05 & prev_result$coefficients3 > 0, ]
prev_result = prev_result[prev_result$V10 == "SKCM", ]

full = get_pathways(prev_result)
gene_list = get_genes(pathways, full)
gene_list = unique(gene_list)

full = strsplit2(pathways, "\t", 1)

p_vals = rep(NA, length(pathways))
coef = rep(NA, length(pathways))
pathway = rep(NA, length(pathways))
num = rep(NA, length(pathways))
clonal_GV = rep(NA, length(pathways))
clonal_no_GV = rep(NA, length(pathways))

temp = MUT2[is.element(merged$V1, gene_list), ]

temp = apply(temp, FUN=not_zero_v2, MARGIN = 2)
#temp = colSums(temp)
sum = sum(temp)
test = as.data.frame(cbind(temp, colnames(MUT2)), stringsAsFactors=FALSE)

meta2 = meta[, 1:2]
colnames(test) = c("GV", "Run")
merged_test = merge(test, meta2, by = "Run")
merged_test2 = merge(merged_test, clinical, by = "Patient")
merged_test2 = merged_test2[!merged_test2$BR == "MR", ]
merged_test2$BR = factor(merged_test2$BR, levels = c("PD", "SD", "PR", "CR"))

require(foreign)
require(ggplot2)
require(MASS)
require(Hmisc)
require(reshape2)
library(survival)

	cox.test = coxph(
		Surv(merged_test2$PFS, merged_test2$progressed) ~ GV + merged_test2[,11] + Tx + Brain_Met + LN_Met + Lung_Met + Liver_Visc_Met + Bone_Met + 
		priorMAPKTx + priorCTLA4 + total_muts, data = merged_test2)

	tab = summary(cox.test)$coefficients

	coxph(
		Surv(merged_test2$PFS, merged_test2$progressed) ~ GV, data = merged_test2)


	library(survminer)
	fit <- survfit(Surv(merged_test2$PFS/365, merged_test2$progressed) ~ merged_test2$GV)
	# Visualize with survminer
	pdf("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Results/Validation/PFS.pdf", height = 5, width = 5)
	ggsurvplot(fit, data = merged_test2, xlab = "Years")
	dev.off()


	merged_test2$GV = as.numeric(merged_test2$GV)

	m <- polr(BR ~ GV + merged_test2[,11] + Tx + Brain_Met + LN_Met + Lung_Met + Liver_Visc_Met + Bone_Met + 
		priorMAPKTx + priorCTLA4, data = merged_test2, Hess = TRUE)

	(ctable <- coef(summary(m)))
	p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
	tab = (ctable <- cbind(ctable, "p value" = p))

	model <- lm(clonal_muts ~ GV + merged_test2[,11] + Tx + Brain_Met + LN_Met + Lung_Met + Liver_Visc_Met + Bone_Met + 
		priorMAPKTx + priorCTLA4, data = merged_test2) 


	m <- polr(BR ~ GV, data = merged_test2, Hess = TRUE)

	(ctable <- coef(summary(m)))
	p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
	tab = (ctable <- cbind(ctable, "p value" = p))

	summary(model)

	wilcox.test(merged_test2$total_muts[merged_test2$GV == 1], merged_test2$total_muts[merged_test2$GV == 0])
	wilcox.test(merged_test2$nonsyn_muts[merged_test2$GV == 1], merged_test2$nonsyn_muts[merged_test2$GV == 0])
	wilcox.test(merged_test2$clonal_muts[merged_test2$GV == 1], merged_test2$clonal_muts[merged_test2$GV == 0])


 # library
library(ggplot2)
tab = table(merged_test2$BR, merged_test2$GV)
frame = as.data.frame(tab, stringsAsFactors=FALSE)
colnames(frame) = c("BR", "PGV", "Pt")
div = c(rep(sum(frame$Pt[frame$PGV == 0]), 4), rep(sum(frame$Pt[frame$PGV == 1]), 4))
frame$Pt = frame$Pt / div * 100

frame$BR = factor(frame$BR, levels = c("CR", "PR", "SD", "PD"))

frame = frame[order(frame$BR), ]

pdf("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Results/Validation/BR.pdf", height = 5, width = 5)
p <- ggplot(frame, aes(fill=BR, y=Pt, x=PGV)) + 
    geom_bar(position="stack", stat="identity") + ylab("Percentage") + theme_bw() + 
    scale_x_discrete(labels=c("No", "Yes")) + xlab("Pathogenic Germline Variant")
print(p)
dev.off()



# POWER ANALYSIS # 

library(powerSurvEpi)
library(survival)

pathways = readLines("/project/dutta_tumor/Ajay/PATH_GV_SIG/Data/Gene_Sets/c2.cp.reactome.v7.0.symbols.gmt")

prev_result = fread3("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Results/Discovery/ind_cancer_pathway.txt", sep = '\t')
prev_result = prev_result[prev_result$p3_adj < 0.05 & prev_result$coefficients3 > 0, ]
prev_result = prev_result[prev_result$V10 == "SKCM", ]
protected = prev_result

pow = rep(NA, nrow(protected))

for(i in 1:nrow(protected)) {
	prev_result = protected
	prev_result = prev_result[i, ]

	full = get_pathways(prev_result)
	gene_list = get_genes(pathways, full)
	gene_list = unique(gene_list)

	temp = MUT2[is.element(merged$V1, gene_list), ]
	temp = apply(temp, FUN=not_zero_v2, MARGIN = 2)
	test = as.data.frame(cbind(temp, colnames(MUT2)), stringsAsFactors=FALSE)

	meta2 = meta[, 1:2]
	colnames(test) = c("GV", "Run")
	merged_test = merge(test, meta2, by = "Run")
	merged_test2 = merge(merged_test, clinical, by = "Patient")
	merged_test2 = merged_test2[!merged_test2$BR == "MR", ]
	merged_test2$BR = factor(merged_test2$BR, levels = c("PD", "SD", "PR", "CR"))

	power_frame = merged_test2[, c(15,31,3)]
	power_frame$group = "E"
	power_frame$group[power_frame$GV == 0] = "C"
	power_frame = power_frame[,c(1,2,4)]
	colnames(power_frame) = c("times", "status", "group")
	power_frame$group = factor(power_frame$group)
	C = sum(power_frame$group == "C")
	E = sum(power_frame$group == "E")
	threshold = 0.05/nrow(protected)

	pow[i] = powerCT(formula = Surv(times, status) ~ group, dat = power_frame, nE = E, nC = C, RR = 2, alpha = threshold)$power

}

# PFS with only disease gene set #

pathways = readLines("/project/dutta_tumor/Ajay/PATH_GV_SIG/Data/Gene_Sets/c2.cp.reactome.v7.0.symbols.gmt")

prev_result = fread3("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Results/Discovery/ind_cancer_pathway.txt", sep = '\t')
prev_result = prev_result[prev_result$p3_adj < 0.05 & prev_result$coefficients3 > 0, ]
prev_result = prev_result[prev_result$V10 == "SKCM", ]
prev_result = prev_result[2, ]

full = get_pathways(prev_result)
gene_list = get_genes(pathways, full)
gene_list = unique(gene_list)

p_vals = rep(NA, length(pathways))
coef = rep(NA, length(pathways))
pathway = rep(NA, length(pathways))
num = rep(NA, length(pathways))
clonal_GV = rep(NA, length(pathways))
clonal_no_GV = rep(NA, length(pathways))

temp = MUT2[is.element(merged$V1, gene_list), ]

temp = apply(temp, FUN=not_zero_v2, MARGIN = 2)
#temp = colSums(temp)
sum = sum(temp)
test = as.data.frame(cbind(temp, colnames(MUT2)), stringsAsFactors=FALSE)

meta2 = meta[, 1:2]
colnames(test) = c("GV", "Run")
merged_test = merge(test, meta2, by = "Run")
merged_test2 = merge(merged_test, clinical, by = "Patient")
merged_test2 = merged_test2[!merged_test2$BR == "MR", ]
merged_test2$BR = factor(merged_test2$BR, levels = c("PD", "SD", "PR", "CR"))

require(foreign)
require(ggplot2)
require(MASS)
require(Hmisc)
require(reshape2)
library(survival)

	cox.test = coxph(
		Surv(merged_test2$PFS, merged_test2$progressed) ~ GV + merged_test2[,11] + Tx + Brain_Met + LN_Met + Lung_Met + Liver_Visc_Met + Bone_Met + 
		priorMAPKTx + priorCTLA4 + total_muts, data = merged_test2)

	tab = summary(cox.test)$coefficients

	coxph(
		Surv(merged_test2$PFS, merged_test2$progressed) ~ GV, data = merged_test2)


	library(survminer)
	fit <- survfit(Surv(merged_test2$PFS/365, merged_test2$progressed) ~ merged_test2$GV)
	# Visualize with survminer
	pdf("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Results/Validation/PFS_disease.pdf", height = 5, width = 5)
	ggsurvplot(fit, data = merged_test2, xlab = "Years")
	dev.off()


	merged_test2$GV = as.numeric(merged_test2$GV)

	m <- polr(BR ~ GV + merged_test2[,11] + Tx + Brain_Met + LN_Met + Lung_Met + Liver_Visc_Met + Bone_Met + 
		priorMAPKTx + priorCTLA4, data = merged_test2, Hess = TRUE)

	(ctable <- coef(summary(m)))
	p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
	tab = (ctable <- cbind(ctable, "p value" = p))

	model <- lm(clonal_muts ~ GV + merged_test2[,11] + Tx + Brain_Met + LN_Met + Lung_Met + Liver_Visc_Met + Bone_Met + 
		priorMAPKTx + priorCTLA4, data = merged_test2) 


	m <- polr(BR ~ GV, data = merged_test2, Hess = TRUE)

	(ctable <- coef(summary(m)))
	p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
	tab = (ctable <- cbind(ctable, "p value" = p))

	summary(model)

	wilcox.test(merged_test2$total_muts[merged_test2$GV == 1], merged_test2$total_muts[merged_test2$GV == 0])
	wilcox.test(merged_test2$nonsyn_muts[merged_test2$GV == 1], merged_test2$nonsyn_muts[merged_test2$GV == 0])
	wilcox.test(merged_test2$clonal_muts[merged_test2$GV == 1], merged_test2$clonal_muts[merged_test2$GV == 0])

library(ggplot2)
tab = table(merged_test2$BR, merged_test2$GV)
frame = as.data.frame(tab, stringsAsFactors=FALSE)
colnames(frame) = c("BR", "PGV", "Pt")
div = c(rep(sum(frame$Pt[frame$PGV == 0]), 4), rep(sum(frame$Pt[frame$PGV == 1]), 4))
frame$Pt = frame$Pt / div * 100

frame$BR = factor(frame$BR, levels = c("CR", "PR", "SD", "PD"))

frame = frame[order(frame$BR), ]

pdf("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Results/Validation/BR_disease.pdf", height = 5, width = 5)
p <- ggplot(frame, aes(fill=BR, y=Pt, x=PGV)) + 
    geom_bar(position="stack", stat="identity") + ylab("Percentage") + theme_bw() + 
    scale_x_discrete(labels=c("No", "Yes")) + xlab("Pathogenic Germline Variant")
print(p)
dev.off()