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

setwd("/project/dutta_tumor/Ajay/Coverage/dbGAP/SRA/dbGaP-24960/fastq/")
DIR = list.files()
main = "/project/dutta_tumor/Ajay/Coverage/dbGAP/SRA/dbGaP-24960/fastq/"
list = list()

for(i in 1:length(DIR)) {
	newdir = paste(main, DIR[i], sep = '')
	file1 = paste(newdir, "/", DIR[i], "_1.fastq", sep = '')
	file2 = paste(newdir, "/", DIR[i], "_1.fastq", sep = '')
	res = as.data.frame(cbind(DIR[i], "lib1", "1", file1, file2), stringsAsFactors=FALSE)
	list[[i]] = res
}

frame = as.data.frame(do.call("rbind", list), stringsAsFactors=FALSE)
library(data.table)
fwrite(frame, "/project/dutta_tumor/Ajay/Coverage/dbGAP/align_cromwell.txt", sep = '\t', col.names = FALSE)


