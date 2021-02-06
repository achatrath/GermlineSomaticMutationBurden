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

p1 = rep(NA, length(genes))
coefficients1 = rep(NA, length(genes))

p2 = rep(NA, length(genes))
coefficients2 = rep(NA, length(genes))

p3 = rep(NA, length(genes))
coefficients3 = rep(NA, length(genes))

num = rep(NA, length(genes))
list = list()

gene_name = rep(NA, length(genes))
cancers = rep(NA, length(genes))

for(i in 1:length(genes)) {
	temp = PGV2[is.element(PGV2$HUGO_Symbol, genes[i]), ]

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
	table = table(temp$cancer)
	cancers = paste(paste(names(table), " (", as.numeric(table), ")", sep = ''), collapse = ", ")
}

result = as.data.frame(cbind(p1, coefficients1, p2, coefficients2, p3, coefficients3, num, cancers, genes), stringsAsFactors=FALSE)
result = result[complete.cases(result), ]
result$p1_adj = p.adjust(result$p1, method = "BH")
result$p2_adj = p.adjust(result$p2, method = "BH")
result$p3_adj = p.adjust(result$p3, method = "BH")
result[result$p3_adj < 0.05, ]
result[result$p3 < 0.05, ]
result = result[order(result$p3_adj), ]

fwrite(result, "/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Results/Discovery/pancancer_gene_result.txt", sep = '\t')

# Make Plot #

# Clonal nonsyn #
result$p = as.numeric(result$p3)
result$val = -log(as.numeric(result$p))
result$Effect = as.numeric(result$coefficients3)
intercept1 = -log(0.05)
intercept2 = -log(0.05/nrow(result))
result$Label = ""
result$Label[result$p < 0.05] = result$genes[result$p < 0.05]

library(ggplot2)
library(ggrepel)
pdf("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Results/Discovery/Clonal_Nonsyn_pancancer_gene.pdf", height = 5, width = 8)
p <- ggplot(result, aes(x=Effect, y=val)) + geom_point() +
   geom_hline(yintercept=intercept1, linetype="dashed", color = "green", size=2) + 
   geom_hline(yintercept=intercept2, linetype="dashed", color = "red", size=2)

p <- p + geom_label_repel(aes(label = Label), box.padding = 0.35, point.padding = 0.5, segment.color = "grey50")
p <- p + theme_bw(base_size=16) + xlab("Additional Clonal Nonsynonymous Somatic Mutations Per MB") + ylab("-log(p)")
p
dev.off()









