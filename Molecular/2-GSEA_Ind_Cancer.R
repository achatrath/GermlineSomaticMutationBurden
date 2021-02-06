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

get_pathways_list <- function(pathways) {
	vec <- rep(NA, length(pathways))
	for(i in 1:length(pathways)) {
		vec[i] = unlist(strsplit(pathways[[i]], "\t"))[1]
	}
	return(vec)
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

# full = full[, medians > 1]
# full2 = scale(full)

rm(PGV)
rm(CDR)
rm(list)
gc()

library(biomaRt)
ensg_id = strsplit2(colnames(full), '[.]', 1)

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
result = getBM(attributes=c("ensembl_gene_id", "external_gene_name"), filters = c("ensembl_gene_id"),
	values = ensg_id, mart = mart)

prev_result = fread3("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Results/Discovery/ind_cancer_pathway.txt", sep = '\t')
prev_result = prev_result[prev_result$p3_adj < 0.05, ]
keep_pathways <- c(1,3,12,14,17,20,22,23,28,29)
prev_result = prev_result[keep_pathways, ]
reformatted_prev = reformat_prev_result(prev_result)

pathways = readLines("/project/dutta_tumor/Ajay/PATH_GV_SIG/Data/Gene_Sets/c2.cp.reactome.v7.0.symbols.gmt")
pathway_list = get_pathways_list(pathways)

# keep_pathways <- prev_result$pathway
# > keep_pathways
#  [1] "REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_DNA_REPAIR_GENES"                                                                                                                                                                                                                                                                                                 
#  [3] "REACTOME_TRANSCRIPTIONAL_REGULATION_BY_TP53"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
# [12] "REACTOME_G2_M_CHECKPOINTS, REACTOME_G2_M_DNA_DAMAGE_CHECKPOINT, REACTOME_REGULATION_OF_TP53_ACTIVITY, REACTOME_REGULATION_OF_TP53_ACTIVITY_THROUGH_PHOSPHORYLATION"                                                                                                                                                             
# [14] "REACTOME_CELL_CYCLE_CHECKPOINTS"                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
# [17] "REACTOME_CELL_CYCLE"                                                                                                                                                                                                                                                                                                                                              
# [20] "REACTOME_CELL_CYCLE"                                                                                                                                                                                                                                                                             
# [22] "REACTOME_CELL_CYCLE_CHECKPOINTS"                                                                                                                                                          
# [23] "REACTOME_REGULATION_OF_TP53_ACTIVITY"                                                                                                                                                                                                                                                            
# [28] "REACTOME_G2_M_CHECKPOINTS, REACTOME_G2_M_DNA_DAMAGE_CHECKPOINT, REACTOME_REGULATION_OF_TP53_ACTIVITY_THROUGH_PHOSPHORYLATION"                                                             
# [29] "REACTOME_TRANSCRIPTIONAL_REGULATION_BY_TP53"                                                                                                                                              

for(i in 1:nrow(reformatted_prev)) {
	temp = merge1
	temp = temp[temp$type == reformatted_prev$cancer[i], ]
	temp$path_variant = 0
	index = which(pathway_list == reformatted_prev$pathway[i])
	genes = unlist(strsplit(pathways[index], "\t"))

	full2 = full[is.element(rownames(full), temp$pt_barcode), ]
	medians = apply(full2, FUN=median, MARGIN=2)
	full2 = full2[, medians > 1]
	full2 = scale(full2)

	temp_PGV2 = PGV2[is.element(PGV2$HUGO_Symbol, genes), ]
	temp$path_variant[is.element(temp$bcr_patient_barcode, temp_PGV2$bcr_patient_barcode)] = 1
	estimates = list()

	library(doParallel)
	library(doSNOW)
	cl <- makeCluster(7, outfile = "")
	registerDoParallel(cl)

	estimates <- foreach(j=1:ncol(full2)) %dopar% {
		a = as.data.frame(cbind(rownames(full2), full2[,j]), stringsAsFactors = FALSE)
		colnames(a) = c("pt_barcode", colnames(full2)[j])
		b = temp

		c = merge(a, b, by = "pt_barcode")

		y = as.numeric(c[,2])

		if(reformatted_prev$cancer[i] == "UCEC")
			model = lm(y ~ c$path_variant + c$age_at_initial_pathologic_diagnosis +
				c[,13] + c[,14] + c[,15] + c[,16])
		else 
			model = lm(y ~ c$path_variant + c$age_at_initial_pathologic_diagnosis +
				c$gender + c[,13] + c[,14] + c[,15] + c[,16])
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
	setwd("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Results/Molecular/GSEA/Ind_Cancer_Pathway")
	file_name = paste(reformatted_prev[i,1], reformatted_prev[i, 2], sep = "-")
	fwrite(res3, file_name, col.names = FALSE, sep = '\t')
}




