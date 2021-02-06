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

setwd("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/TCGA_FILES")

PGV = fread3("PCA_pathVar_integrated_filtered_adjusted.tsv", sep = '\t')
CDR = fread3('TCGA-CDR.csv', sep = '\t')
TCGAA = fread3('STRUCTURE.txt', sep = '\t')


PGV2 = PGV[,c(1:3)]
PGV2$pt_barcode = strsplit2(PGV2[,1], "-", 3)
PGV2$path_variant = 1

CDR2 = CDR[,c(2:5)]
CDR2$pt_barcode = strsplit2(CDR2[,1], "-", 3)
TCGAA$pt_barcode = strsplit2(TCGAA[,1], "-", 3)

merge1 = merge(CDR2, TCGAA, by = "pt_barcode")

setwd("/project/dutta_tumor/Ajay/BACKUP/ac4ec/CGC/Data/RNA-Seq-Compiled/Tumor")
FILES = list.files()
list = list()

for(i in 1:length(FILES)) {
	print(i)
	data = t(fread2(FILES[i], sep = '\t'))
	list[[i]] = data
	pt = c(pt, rep(FILES[i], nrow(data)))
}

full = as.data.frame(do.call("rbind", list), stringsAsFactors=FALSE)
medians = apply(full, FUN=median, MARGIN=2)

full = full[, medians > 1]
full2 = scale(full)

#quantile(apply(full2, FUN=mean, MARGIN=2))
#quantile(apply(full2, FUN=sd, MARGIN=2))

#quantile(apply(full2, FUN=mean, MARGIN=1))
#quantile(apply(full2, FUN=sd, MARGIN=1))

rm(PGV)
rm(CDR)
gc()

library(biomaRt)
ensg_id = strsplit2(colnames(full2), '[.]', 1)

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
result = getBM(attributes=c("ensembl_gene_id", "external_gene_name"), filters = c("ensembl_gene_id"),
	values = ensg_id, mart = mart)

# genes = c("APC", "ERCC3", "FANCL", "MSH2", "PMS2", "SLC25A13", "TP53")
genes = c("APC")

for(i in 1:length(genes)) {
	temp = merge1
	temp$path_variant = 0
	temp_PGV2 = PGV2[PGV2$HUGO_Symbol == genes[i], ]
	temp$path_variant[is.element(temp$bcr_patient_barcode, temp_PGV2$bcr_patient_barcode)] = 1
	estimates = list()

	library(doParallel)
	library(doSNOW)
	cl <- makeCluster(20, outfile = "")
	registerDoParallel(cl)

	estimates <- foreach(j=1:ncol(full2)) %dopar% {
		a = as.data.frame(cbind(rownames(full2), full2[,j]), stringsAsFactors = FALSE)
		colnames(a) = c("pt_barcode", colnames(full2)[j])
		b = temp

		c = merge(a, b, by = "pt_barcode")

		y = as.numeric(c[,2])

		model = lm(y ~ c$path_variant + c$type + c$age_at_initial_pathologic_diagnosis +
			c$gender + c[,13] + c[,14] + c[,15] + c[,16])
		# model = lm(y ~ c$path_variant + c$type)
		estimate = summary(model)$coefficients[2,1]
		p_val = summary(model)$coefficients[2,4]

		temp_res = as.data.frame(cbind(estimate, p_val), stringsAsFactors=FALSE)
		return(temp_res)
	}

	stopCluster(cl) 

	print("NEXT")

	res = as.data.frame(do.call("rbind", estimates), stringsAsFactors=FALSE)

	res = as.data.frame(cbind(colnames(full2), res), stringsAsFactors=FALSE)
	res = res[order(as.numeric(res$estimate)), ]
	res[,1] = as.character(res[,1])
	res$ensembl_gene_id = strsplit2(res[,1], "[.]", 1)
	res2 = merge(res, result, by = "ensembl_gene_id")
	options(scipen=999)
	res3 = res2[, c(5, 3)]
	res3$estimate = as.numeric(res3$estimate)
	res3 = res3[order(res3$estimate, decreasing = TRUE), ]
	setwd("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Results/Molecular/GSEA/Gene")
	file_name = paste(genes[i], ".rnk", sep = '')
	fwrite(res3, file_name, col.names = FALSE, sep = '\t')
}