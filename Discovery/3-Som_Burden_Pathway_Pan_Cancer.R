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

parse_contrib <- function(contr, CDR) {
  pt = strsplit2(unlist(strsplit(contr, ", ")), "_", 1)
  cancers = sort(unique(CDR$type))
  vector = rep(NA, length(cancers))

  for(i in 1:length(cancers)) {
    temp = CDR[CDR$type == cancers[i], ]
    vector[i] = sum(is.element(pt, temp$bcr_patient_barcode))
  }

  vector = as.data.frame(t(as.data.frame(vector, stringsAsFactors=FALSE)), stringsAsFactors=FALSE)
  colnames(vector) = cancers
  return(vector)
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

tab = table(PGV2$HUGO_Symbol)
tab2 = tab[tab >= 5]

genes = names(tab[tab >= 5])

p1 = rep(NA, length(pathways))
coefficients1 = rep(NA, length(pathways))

p2 = rep(NA, length(pathways))
coefficients2 = rep(NA, length(pathways))

p3 = rep(NA, length(pathways))
coefficients3 = rep(NA, length(pathways))

num = rep(NA, length(pathways))
list = list()
pathway = rep(NA, length(pathways))
cancer_list = rep(NA, length(pathways))

pathways = readLines("/project/dutta_tumor/Ajay/PATH_GV_SIG/Data/Gene_Sets/c2.cp.reactome.v7.0.symbols.gmt")

perturbed_genes = rep(NA, length(pathways))

for(i in 1:length(pathways)) {
	genes = unlist(strsplit(pathways[[i]], '\t'))
	temp = PGV2[is.element(PGV2$HUGO_Symbol, genes), ]

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

	tab = sort(table(temp$cancer))
	cancer_list[i] = paste(paste(names(tab), " (", tab, ")", sep = ''), collapse = ', ')
}

result = as.data.frame(cbind(p1, coefficients1, p2, coefficients2, p3, coefficients3, num, pathway, perturbed_genes, cancer_list), stringsAsFactors=FALSE)
result = result[complete.cases(result), ]
result = overlapped_pathways(result)
result$p1_adj = p.adjust(result$p1, method = "BH")
result$p2_adj = p.adjust(result$p2, method = "BH")
result$p3_adj = p.adjust(result$p3, method = "BH")

result$gene_list = unlist(lapply(result$perturbed_genes, FUN=get_genes))

fwrite(result, "/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Results/Discovery/pancancer_pathway.txt", sep = '\t')

# Make Plot #

colors = fread3("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/TCGA_FILES/colors.tsv", sep = '\t')
col = unlist(strsplit(colors$V1[33], ";"))

library(parallel)
data = result
effect = data$coefficients3
h = -log(data$p3_adj, 10)
contr = mclapply(data$perturbed_genes, CDR = CDR, FUN=parse_contrib, mc.cores = 7)
contr = as.data.frame(do.call("rbind", contr), stringsAsFactors=FALSE)

frame = as.data.frame(cbind(effect, h, contr), stringsAsFactors=FALSE)
frame$effect = as.numeric(frame$effect)
frame$scaled_effect = as.numeric(scale(frame$effect))

library(circlize)
library(gtools)
library(ggplot2)
library(survival)
library(survminer)
library(scatterpie)

frame$region = factor(1:nrow(frame))

setwd("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Results/Discovery")
pdf("Scatterpie.pdf", height = 10, width = 15)
ggplot() + geom_scatterpie(aes(x=scaled_effect, y=h, group=region, r = 0.3), data=frame,
                           cols=colnames(frame)[3:35]) + coord_equal() + 
	scale_fill_manual(values=col) + theme_bw(base_size=18) + geom_hline(yintercept = -log(0.05, 10), lwd = 1.2) + xlab("Scaled Effect Size") + 
	ylab("")+ labs(fill = "Cancer")
dev.off()



