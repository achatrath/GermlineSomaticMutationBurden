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

medians = apply(full, FUN=median, MARGIN=2)
full = full[, medians > 1]
full2 = scale(full)

rm(PGV)
rm(CDR)
rm(list)
gc()
rm(full)

library(biomaRt)
ensg_id = strsplit2(colnames(full2), '[.]', 1)

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
result = getBM(attributes=c("ensembl_gene_id", "external_gene_name"), filters = c("ensembl_gene_id"),
	values = ensg_id, mart = mart)

prev_result = fread3("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Results/Discovery/pancancer_pathway.txt", sep = '\t')
prev_result = prev_result[prev_result$p3_adj < 0.05, ]
keep_pathways <- c(1,2,4,5,6,7,8)
prev_result = prev_result[keep_pathways, ]
reformatted_prev = reformat_prev_result(prev_result)

pathways = readLines("/project/dutta_tumor/Ajay/PATH_GV_SIG/Data/Gene_Sets/c2.cp.reactome.v7.0.symbols.gmt")
pathway_list = get_pathways_list(pathways)

# keep_pathways <- prev_result$pathway
# # > keep_pathways
#  [1] "REACTOME_DEGRADATION_OF_BETA_CATENIN_BY_THE_DESTRUCTION_COMPLEX"                                                                                                                                                                                                                           
#  [2] "REACTOME_BETA_CATENIN_PHOSPHORYLATION_CASCADE, REACTOME_DISASSEMBLY_OF_THE_DESTRUCTION_COMPLEX_AND_RECRUITMENT_OF_AXIN_TO_THE_MEMBRANE, REACTOME_SIGNALING_BY_WNT_IN_CANCER, REACTOME_PHOSPHORYLATION_SITE_MUTANTS_OF_CTNNB1_ARE_NOT_TARGETED_TO_THE_PROTEASOME_BY_THE_DESTRUCTION_COMPLEX"                                                                                                                                                                                                                                             
#  [4] "REACTOME_DEACTIVATION_OF_THE_BETA_CATENIN_TRANSACTIVATING_COMPLEX"                                                                                                                                                                                                                         
#  [5] "REACTOME_PROGRAMMED_CELL_DEATH"                                                                                                                                                                                                                                                            
#  [6] "REACTOME_REGULATION_OF_KIT_SIGNALING"                                                                                                                                                                                                                                                      
#  [7] "REACTOME_APOPTOTIC_CLEAVAGE_OF_CELLULAR_PROTEINS, REACTOME_APOPTOTIC_EXECUTION_PHASE"                                                                                                                                                                                                      
#  [8] "REACTOME_SIGNALING_BY_WNT, REACTOME_TCF_DEPENDENT_SIGNALING_IN_RESPONSE_TO_WNT"                                                                                                                                                                                                            
                                                                                                                                                

for(i in 1:nrow(reformatted_prev)) {
	temp = merge1
	temp$path_variant = 0
	index = which(pathway_list == reformatted_prev$pathway[i])
	genes = unlist(strsplit(pathways[index], "\t"))

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

		model = lm(y ~ c$path_variant + c$age_at_initial_pathologic_diagnosis + c$type + 
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
	setwd("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Results/Molecular/GSEA/Pancancer_Pathway")
	file_name = paste(reformatted_prev[i,1], ".rnk", sep = "")
	fwrite(res3, file_name, col.names = FALSE, sep = '\t')
}




