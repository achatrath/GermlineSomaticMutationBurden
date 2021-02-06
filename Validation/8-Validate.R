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
		list[[i]] = unlist(strsplit(prev_result$V2[i], ", "))
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


# individual genes not testable #
pathways = readLines("/project/dutta_tumor/Ajay/PATH_GV_SIG/Data/Gene_Sets/c2.cp.reactome.v7.0.symbols.gmt")

p_vals = rep(NA, length(pathways))
coef = rep(NA, length(pathways))
pathway = rep(NA, length(pathways))
num = rep(NA, length(pathways))

for(i in 1:length(pathways)) {
	split = unlist(strsplit(pathways[i], "\t"))
	temp = MUT2[is.element(merged$V1, split), ]

	if(class(temp) == "matrix") {
		if(nrow(temp) == 0) {
			next
		}
	}


	if(class(temp) == 'numeric') {
		temp = as.numeric(temp != 0)
		sum = sum(temp)
		test = as.data.frame(cbind(temp, colnames(MUT2)), stringsAsFactors=FALSE)
		if(sum < 5) {
			next
		}
	}

	if(class(temp) == "matrix") {
		print(i)
		temp = apply(temp, FUN=not_zero_v2, MARGIN = 2)
		sum = sum(temp)
		test = as.data.frame(cbind(temp, colnames(MUT2)), stringsAsFactors=FALSE)
		if(sum<5) 
			next
	}

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

	#m <- polr(BR ~ GV + merged_test2[,11] + Tx + Brain_Met + LN_Met + Lung_Met + Liver_Visc_Met + Bone_Met + postCTLA4, data = merged_test2, Hess=TRUE)
	m <- polr(BR ~ GV + merged_test2[,11] + merged_test2$Tx + priorMAPKTx + priorCTLA4 + postCTLA4 + postMAPKTx + postCombinedCTLA_PD1, data = merged_test2, Hess=TRUE)
	(ctable <- coef(summary(m)))
	p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
	tab = (ctable <- cbind(ctable, "p value" = p))

	p_vals[i] = as.numeric(tab[1,4])
	coef[i] = as.numeric(exp(tab[1,1]))
	pathway[i] = split[1]
	num[i] = as.numeric(sum(temp))
}

result_frame = as.data.frame(cbind(p_vals, coef, pathway, num), stringsAsFactors=FALSE)
result_frame = result_frame[complete.cases(result_frame), ]
result_frame$p_vals = as.numeric(result_frame$p_vals)
result_frame = result_frame[order(result_frame$p_vals), ]
result_frame = result_frame[result_frame$num >= 5, ]

fwrite(result_frame, "/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Result_v4.txt", sep = '\t')

# Try with PFI #

# individual genes not testable #
pathways = readLines("/project/dutta_tumor/Ajay/PATH_GV_SIG/Data/Gene_Sets/c2.cp.reactome.v7.0.symbols.gmt")

p_vals = rep(NA, length(pathways))
coef = rep(NA, length(pathways))
pathway = rep(NA, length(pathways))
num = rep(NA, length(pathways))
clonal_GV = rep(NA, length(pathways))
clonal_no_GV = rep(NA, length(pathways))

for(i in 1:length(pathways)) {
	split = unlist(strsplit(pathways[i], "\t"))
	temp = MUT2[is.element(merged$V1, split), ]

	if(class(temp) == "matrix") {
		if(nrow(temp) == 0) {
			next
		}
	}


	if(class(temp) == 'numeric') {
		temp = as.numeric(temp != 0)
		sum = sum(temp)
		test = as.data.frame(cbind(temp, colnames(MUT2)), stringsAsFactors=FALSE)
		if(sum < 5) {
			next
		}
	}

	if(class(temp) == "matrix") {
		print(i)
		temp = apply(temp, FUN=not_zero_v2, MARGIN = 2)
		sum = sum(temp)
		test = as.data.frame(cbind(temp, colnames(MUT2)), stringsAsFactors=FALSE)
		if(sum<5) 
			next
	}

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
		Surv(merged_test2$PFS, merged_test2$progressed) ~ GV + merged_test2[,11] + Tx + Tx + Brain_Met + LN_Met + Lung_Met + Liver_Visc_Met + Bone_Met + 
		priorMAPKTx + priorCTLA4 + postCTLA4 + postMAPKTx + postCombinedCTLA_PD1,, data = merged_test2)

	tab = summary(cox.test)$coefficients

	p_vals[i] = as.numeric(tab[1,5])
	coef[i] = as.numeric(tab[1,2])
	pathway[i] = split[1]
	num[i] = as.numeric(sum(temp))
	clonal_GV[i] = mean(merged_test2$clonal_muts[merged_test2$GV == 1])
	clonal_no_GV[i] = mean(merged_test2$clonal_muts[merged_test2$GV == 0])

}

result_frame = as.data.frame(cbind(p_vals, coef, pathway, num, clonal_GV, clonal_no_GV), stringsAsFactors=FALSE)
result_frame = result_frame[complete.cases(result_frame), ]
result_frame$p_vals = as.numeric(result_frame$p_vals)
result_frame = result_frame[order(result_frame$p_vals), ]
result_frame = result_frame[result_frame$num >= 5, ]

fwrite(result_frame, "/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Result_PFS.txt", sep = '\t')

# Try with PFI but with summing #

# Try with PFI #

# individual genes not testable #
pathways = readLines("/project/dutta_tumor/Ajay/PATH_GV_SIG/Data/Gene_Sets/c2.cp.reactome.v7.0.symbols.gmt")

p_vals = rep(NA, length(pathways))
coef = rep(NA, length(pathways))
pathway = rep(NA, length(pathways))
num = rep(NA, length(pathways))
clonal_GV = rep(NA, length(pathways))
clonal_no_GV = rep(NA, length(pathways))

for(i in 1:length(pathways)) {
	split = unlist(strsplit(pathways[i], "\t"))
	temp = MUT2[is.element(merged$V1, split), ]

	if(class(temp) == "matrix") {
		if(nrow(temp) == 0) {
			next
		}
	}


	if(class(temp) == 'numeric') {
		sum = sum(temp != 0)
		test = as.data.frame(cbind(temp, colnames(MUT2)), stringsAsFactors=FALSE)
		if(sum < 5) {
			next
		}
	}

	if(class(temp) == "matrix") {
		print(i)
		temp = colSums(temp)
		sum = sum(temp != 0)
		test = as.data.frame(cbind(temp, colnames(MUT2)), stringsAsFactors=FALSE)
		if(sum<5) 
			next
	}

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
		Surv(merged_test2$PFS, merged_test2$progressed) ~ GV + merged_test2[,11] + Tx + Tx + Brain_Met + LN_Met + Lung_Met + Liver_Visc_Met + Bone_Met + 
		priorMAPKTx + priorCTLA4 + postCTLA4 + postMAPKTx + postCombinedCTLA_PD1,, data = merged_test2)

	tab = summary(cox.test)$coefficients

	p_vals[i] = as.numeric(tab[1,5])
	coef[i] = as.numeric(tab[1,2])
	pathway[i] = split[1]
	num[i] = as.numeric(sum(temp))
	clonal_GV[i] = mean(merged_test2$clonal_muts[merged_test2$GV == 1])
	clonal_no_GV[i] = mean(merged_test2$clonal_muts[merged_test2$GV == 0])

}

result_frame = as.data.frame(cbind(p_vals, coef, pathway, num, clonal_GV, clonal_no_GV), stringsAsFactors=FALSE)
result_frame = result_frame[complete.cases(result_frame), ]
result_frame$p_vals = as.numeric(result_frame$p_vals)
result_frame = result_frame[order(result_frame$p_vals), ]
result_frame = result_frame[result_frame$num >= 5, ]

fwrite(result_frame, "/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Result_PFS.txt", sep = '\t')


# try with PFI but pooling pathways together #
# Try with PFI #

# individual genes not testable #
pathways = readLines("/project/dutta_tumor/Ajay/PATH_GV_SIG/Data/Gene_Sets/c2.cp.reactome.v7.0.symbols.gmt")

prev_result = fread3("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Pan_cancer_pathway_V2.tsv", sep = '\t')
prev_result = prev_result[prev_result$padj < 0.05, ]

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
		priorMAPKTx + priorCTLA4, data = merged_test2)



	tab = summary(cox.test)$coefficients

	merged_test2$GV = as.numeric(merged_test2$GV)

	m <- polr(BR ~ GV + merged_test2[,11] + Tx + Brain_Met + LN_Met + Lung_Met + Liver_Visc_Met + Bone_Met + 
		priorMAPKTx + priorCTLA4, data = merged_test2, Hess = TRUE)
	(ctable <- coef(summary(m)))
	p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
	tab = (ctable <- cbind(ctable, "p value" = p))

	p_vals[i] = as.numeric(tab[1,5])
	coef[i] = as.numeric(tab[1,2])
	pathway[i] = split[1]
	num[i] = as.numeric(sum(temp))
	clonal_GV[i] = mean(merged_test2$clonal_muts[merged_test2$GV == 1])
	clonal_no_GV[i] = mean(merged_test2$clonal_muts[merged_test2$GV == 0])

	t.test(merged_test2$total_muts[merged_test2$GV == 1], merged_test2$total_muts[merged_test2$GV == 0])
	t.test(merged_test2$nonsyn_muts[merged_test2$GV == 1], merged_test2$nonsyn_muts[merged_test2$GV == 0])
	t.test(merged_test2$clonal_muts[merged_test2$GV == 1], merged_test2$clonal_muts[merged_test2$GV == 0])



result_frame = as.data.frame(cbind(p_vals, coef, pathway, num, clonal_GV, clonal_no_GV), stringsAsFactors=FALSE)
result_frame = result_frame[complete.cases(result_frame), ]
result_frame$p_vals = as.numeric(result_frame$p_vals)
result_frame = result_frame[order(result_frame$p_vals), ]
result_frame = result_frame[result_frame$num >= 5, ]

fwrite(result_frame, "/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Result_PFS.txt", sep = '\t')







