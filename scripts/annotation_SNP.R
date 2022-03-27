library("ggplot2")
library("Rsamtools")
library("GenomicAlignments")
library("rtracklayer")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("GenomicRanges")
library("AnnotationHub")
library("knitr")
library("gtools")
library("data.table")
library("stringi")
library(GBJ)
library(metap)
library(multtest)
library(Hmisc)
library(devtools)
library(SNPRelate)
library(gdsfmt)
library(dplyr)

###Loading database
##Download SNP database
#SNP Database 
#http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-52/plants/gff3/sorghum_bicolor/
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Data for sorghum/sorghum/SNP annotation/Sorghum")
snp.db <- read.table(file ="Sbicolor_454_v3.1.1.gene_exons.gff3", sep = "\t", header = FALSE)
colnames(snp.db) <- c("Chromosome","Database","Region","Start","End","NA","Strand","NA2","Gene")

###Loading query SNPs for GWAS from RDS file
##Location in the server: /rsstu/users/r/rrellan/sara/SorghumGEA/results/GLM_20220222 - GLM
##Location in the server: /rsstu/users/r/rrellan/sara/SorghumGEA/results/GLM_20220224 - MLM
#setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Data for sorghum/sorghum/GWAS.results")

for(i in 1:10){
  assign(paste0("query.snp.gwas", i) , readRDS(file = paste0("glm_sol_VL_",i,"_20220224_03_09.RDS")))
}
for(i in 1:10){
  assign(paste0("query.snp.gwas",i), as.data.frame(paste0("query.snp.gwas",i,"$GLM_Stats")))
  assign(paste0("query.snp.gwas",i), paste0("query.snp.gwas",i,"$GLM_Stats")[,c(2,3,4,5,6)])
}

###Subsetting chr from db
for(i in 1:10){
  assign(paste0("db.",i), snp.db[which(snp.db$Chromosome == paste0("Chr",i))])
}


##Making GRanges for Database 
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Data for sorghum/sorghum/GWAS.results")
for(i in 01:10){
  assign(paste0("gr.db", i) , GRanges(seqnames = paste0("chr",i), ranges = IRanges(start = get(paste0("db.",i))[,"Start"], end = get(paste0("db.",i))[,"End"]), strand = get(paste0("db.",i))[,"Strand"], Region = get(paste0("db.",i))[,"Region"], Gene = get(paste0("db.",i))[,"Gene"]))
}

#Making GRanges for gwas Query
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Data for sorghum/sorghum/GWAS.results")
for(i in 1:10){
  assign(paste0("query.snp.gwas",i), readRDS(file = paste0("glm_sol_VL_",i,"_20220222_13_55.RDS")))
}
for(i in 1:10){
  assign(paste0("query.snp.gwas",i), as.data.frame(paste0("query.snp.gwas",i,"$GLM_Stats")))
  assign(paste0("query.snp.gwas",i), paste0("query.snp.gwas",i[,c(2,3,4,5,6)]))
}

for(i in 1:10){
  assign(paste0("gr.q", i) , GRanges(seqnames = paste0("chr",i), ranges = IRanges(start = get(paste0("query",i,".gwas"))[,"Pos"], width =1, fstat = get(paste0("query",i,".gwas"))[,"p"], maker = get(paste0("query",i,".gwas"))[,"Marker"])))
}


##Overlaps
for(i in 1:10){
  assign(paste0("common",i), as.data.frame(findOverlapPairs(paste0("gr.db",i), paste0("gr.q.",i))))
}

