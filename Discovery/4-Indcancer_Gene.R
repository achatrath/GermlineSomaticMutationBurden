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

list = list()
cancers = sort(unique(CDR2$type))

for(i in 1:length(cancers)) {
	temp = PGV2[PGV2$cancer == cancers[i], ]
	tab = table(temp$HUGO_Symbol)
	tab2 = tab[tab >= 5]
	if(length(tab2) == 0)
		next

	genes = names(tab[tab >= 5])
	p3 = rep(NA, length(genes))
	coefficients3 = rep(NA, length(genes))
	num = rep(NA, length(genes))

	for(j in 1:length(genes)) {
		test = merge
		test = test[test$type == cancers[i], ]
		temp2 = temp[temp$HUGO_Symbol == genes[j], ]
		test$GV = 0
		test$GV[is.element(test$patient, temp2$pt_barcode)] = 1
		test = test[complete.cases(test), ]
		if(sum(test$GV) < 5)
			next
		if(cancers[i] == "OV" | cancers[i] == "PRAD" | cancers[i] == "UCEC")
			linearMod <- try(lm(Overall_Per_MB ~ GV + age_at_initial_pathologic_diagnosis +
				European + test[,9] + test[,10] + test[,11], data=test))
		else
			linearMod <- try(lm(Overall_Per_MB ~ GV + age_at_initial_pathologic_diagnosis + gender +
				European + test[,9] + test[,10] + test[,11], data=test))

		num[j] = sum(test$GV == 1)
		p3[j] = summary(linearMod)$coefficients[2,4]
		coefficients3[j] = summary(linearMod)$coefficients[2,1]
	}

	result = as.data.frame(cbind(num, p3, coefficients3, genes), stringsAsFactors=FALSE)
	result = result[complete.cases(result), ]
	result$p3_adj = p.adjust(result$p3, method = "BH")
	list[[i]] = result
}

full = as.data.frame(do.call("rbind", list), stringsAsFactors=FALSE)



	
