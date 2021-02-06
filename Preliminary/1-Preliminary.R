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

getMedian <- function(plot) {
	cancer_list <- unique(plot$Cancer)
	med = rep(NA, length(cancer_list))
	for(i in 1:length(cancer_list)) {
		temp = plot[plot$Cancer == cancer_list[i] & plot$variable == "Overall", ]
		med[i] = median(temp$value)
	}
	res = as.data.frame(cbind(cancer_list, med), stringsAsFactors=FALSE)
	res$med = as.numeric(res$med)
	res = res[order(res$med, decreasing = TRUE), ]
	return(res)
}

TMB = fread3("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Somatic_Burden_Norm_Cov.tsv", sep = '\t')
CDR = fread3("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/TCGA_FILES/TCGA-CDR.csv", sep = '\t')
CDR$Patient = strsplit2(CDR$bcr_patient_barcode, "-", 3)
CDR = CDR[, c(3,35)]
TMB$Patient = strsplit2(TMB$Barcode, "-", 3)

merged = merge(TMB, CDR, by = "Patient")
merged = merged[complete.cases(merged), ]

colors = fread3("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/TCGA_FILES/colors.tsv", sep = '\t')
col = unlist(strsplit(colors$V1[33], ";"))

library(reshape2)
plot = melt(merged, id.vars = c("Patient", "Barcode", "Coverage", "type"), 
	measure.vars = c("Overall_Count", "Nonsyn_Count", "Nonsyn_Clonal_Count", "Overall_Per_MB", "Nonsyn_Per_MB", "Clonal_Nonsyn_Per_MB"))

plot$value2 = log(plot$value)
plot=plot[plot$value2 != -Inf, ]

plot$variable = factor(plot$variable)
levels(plot$variable)
levels(plot$variable) <- c("Overall", "Nonsynonymous", "Clonal Nonsynonymous", 
	"Overall Per MB", "Nonsynonymous Per MB", "Clonal Nonsynonymous Per MB")

colnames(plot)[4] = "Cancer"
med = getMedian(plot)
plot$Cancer = factor(plot$Cancer, levels = med$cancer_list)

library(ggplot2)
pdf("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Results/Preliminary/Clonal_Nonsyn_Per_MB.pdf", height = 20, width = 40)
p <- ggplot(plot, aes(x=Cancer, y=value2, fill = Cancer)) + geom_boxplot(outlier.shape = NA) + facet_wrap(~variable, scales = "free") + 
	theme_bw(base_size = 40) + scale_fill_manual(values = col) + ylab("log(TMB)") + theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
print(p)
dev.off()

library(ggplot2)
plot.f = plot
exclude = c("Overall", "Nonsynonymous", "Clonal Nonsynonymous")
plot.f = plot.f[!is.element(plot.f$variable, exclude), ]
pdf("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Results/Preliminary/Clonal_Nonsyn_Per_MB_filtered.pdf", height = 10, width = 40)
p <- ggplot(plot.f, aes(x=Cancer, y=value2, fill = Cancer)) + geom_boxplot(outlier.shape = NA) + facet_wrap(~variable, scales = "free") + 
	theme_bw(base_size = 40) + scale_fill_manual(values = col) + ylab("log(TMB)") + theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
print(p)
dev.off()

# plot correlogram #
library(RColorBrewer)
mat = merged[,c(3:5,7:9)]
colnames(mat) = c("Overall", "Nonsyn", "Clonal Nonsyn", "Overall/MB", "Nonsyn/MB", "Clonal Nonsyn/MB")
res <- cor(mat, method = "spearman")
colors = brewer.pal(n = 30, name = "RdYlBu")

library(corrplot)
pdf("/project/dutta_tumor/Ajay/Coverage/dbGAP/CHARGER/Results/Preliminary/correlogram.pdf", height = 6, width = 6)
p <- corrplot.mixed(res, tl.pos = "lt", upper = "ellipse", diag = "u")
print(p)
dev.off()










