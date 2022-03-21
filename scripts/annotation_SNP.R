###Packages
library("ggplot2")
library("GenomicRanges")
library("Rsamtools")
library("GenomicAlignments")
library("rtracklayer")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("AnnotationHub")
library("knitr")
library("gtools")
library("data.table")

###Loading database
##download SNP database
#SNP Database 
#http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-52/plants/gff3/sorghum_bicolor/
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Data for sorghum/sorghum/SNP annotation/Sorghum")
snp.db <- read.table(file ="Sbicolor_454_v3.1.1.gene_exons.gff3", sep = "\t", header = FALSE)
colnames(snp.db) <- c("Chromosome","Database","Region","Start","End","NA","Strand","NA2","Gene")


###Loading query SNPs for GWAS from RDS file
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Data for sorghum/sorghum/GWAS.results")
query.snp.gwas1 <- readRDS(file = "glm_sol_VL_01_20220222_13_55.RDS")
query.snp.gwas1 <- as.data.frame(query.snp.gwas1$GLM_Stats)
query.snp.gwas1 <- query.snp.gwas1[,c(2,3,4,6)]


###Subsetting chr from db
db.1 <- snp.db[which(snp.db$Chromosome == "Chr01"), ]
#db.2 <- snp.db[which(snp.db$Chromosome == "chr2"), ]
##db.2 <- db.2[-1,]
#db.3 <- snp.db[which(snp.db$Chromosome == "chr3"), ]
#db.3 <- db.3[-1,]
#db.4 <- snp.db[which(snp.db$Chromosome == "chr4"), ]
#db.4 <- db.4[-1,]
#db.5 <- snp.db[which(snp.db$Chromosome == "chr5"), ]
#db.5 <- db.5[-1,]
#db.6 <- snp.db[which(snp.db$Chromosome == "chr6"), ]
#db.6 <- db.6[-1,]
#db.7 <- snp.db[which(snp.db$Chromosome == "chr7"), ]
#db.7 <- db.7[-1,]
#db.8 <- snp.db[which(snp.db$Chromosome == "chr8"), ]
#db.8 <- db.8[-1,]
#db.9 <- snp.db[which(snp.db$Chromosome == "chr9"), ]
#db.9 <- db.9[-1,]
#db.10 <- snp.db[which(snp.db$Chromosome == "chr10"), ]
#db.10 <- db.10[-1,]



###SNP annotation (Using packages)
#https://kasperdanielhansen.github.io/genbioconductor/html/GenomicRanges_GRanges_Usage.html
#https://combine-australia.github.io/2017-05-19-bioconductor-melbourne/AdvGRanges_Rtracklayer_Rsamtools.html



#Making GRanges for Database 
gr.db1 <- GRanges(seqnames = "chr1", ranges = IRanges(start = db.1$Start, end = db.1$End), strand = db.1$Strand, Region = db.1$Region, Gene = db.1$Gene)
#for(i in 01:10){
#  assign(paste0("gr.db", i) , GRanges(seqnames = paste0("chr",i), ranges = IRanges(start = get(paste0("db.",i))[,"Start"], end = get(paste0("db.",i))[,"End"]), strand = get(paste0("db.",i))[,"Strand"], Region = get(paste0("db.",i))[,"Region"], Gene = get(paste0("db.",i))[,"Gene"]))
#}

#Making GRanges for gwas Query
gr.q.1 <- GRanges(seqnames = "chr1", ranges = IRanges(start = query.snp.gwas1$Pos, width = 1), pvalue = query.snp.gwas1$p, Marker = query.snp.gwas1$Marker)
#for(i in 1:10){
#  assign(paste0("gr.q", i) , GRanges(seqnames = paste0("chr",i), ranges = IRanges(start = get(paste0("query",i,".gwas"))[,"Pos"], width =1, pvalue = get(paste0("query",i,".gwas"))[,"p"], maker = get(paste0("query",i,".gwas"))[,"Marker"])))
#}


#Overlaps
commmon1 <- as.data.frame(findOverlapPairs(gr.db1,gr.q.1))
#commmon2 <- as.data.frame(findOverlapPairs(gr.db2,gr.q2))
#commmon3 <- as.data.frame(findOverlapPairs(gr.db3,gr.q3))
#commmon4 <- as.data.frame(findOverlapPairs(gr.db4,gr.q4))
#commmon5 <- as.data.frame(findOverlapPairs(gr.db5,gr.q5))
#commmon6 <- as.data.frame(findOverlapPairs(gr.db6,gr.q6))
#commmon7 <- as.data.frame(findOverlapPairs(gr.db7,gr.q7))
#commmon8 <- as.data.frame(findOverlapPairs(gr.db8,gr.q8))
#commmon9 <- as.data.frame(findOverlapPairs(gr.db9,gr.q9))
#commmon10 <- as.data.frame(findOverlapPairs(gr.db10,gr.q10))


#Combining pvalues:

#Filtering out only the gene and pvalue columns
gwas <- commmon1
gwas <- gwas[,c(6,7,15,16)]
gwas <- gwas[which(gwas$first.X.Region == "gene"),]
gwas <- gwas[,c(2,3,4)]
colnames(gwas) <- c("Gene","pvalue","Marker")
#Sorting a/c gene name
gwas <- gwas[gtools::mixedorder(gwas$Gene),]


#Table with pvalues
gwas.pvalues <- dcast(setDT(gwas), Gene~rowid(Gene, prefix = "pvalue"), value.var = "pvalue")
gwas.pvalues <- as.data.frame(t(gwas.pvalues))
#Genenames
gene.names <- gwas.pvalues[1,]
#Function
split.names <- function(x,split){
  split.genename <- unlist(strsplit(x, split = ';', fixed = TRUE))[2]
  split.genename2 <- unlist(strsplit(split.genename, split = "=", fixed = TRUE))[2]
  return(split.genename2)
}
gene.names <- apply(gene.names,2,split.names)
gene.names <- as.data.frame(t(gene.names))
#Continue with table with pvalues
gwas.pvalues <- gwas.pvalues[-1,]
gwas.pvalues <- apply(gwas.pvalues,2,sort)
gwas.pvalues <- na.omit(gwas.pvalues)
gwas.pvalues <- stri_list2matrix(gwas.pvalues, byrow = FALSE)
gwas.pvalues <- as.data.frame(gwas.pvalues)
colnames(gwas.pvalues) <- gene.names
#Table with SNP markers
gwas.markers <- dcast(setDT(gwas), Gene~rowid(Gene, prefix = "Marker"), value.var = "Marker")
gwas.markers <- as.data.frame(t(gwas.markers))
gwas.markers <- gwas.markers[-1,]
gwas.markers <- apply(gwas.markers,2,sort)
gwas.markers <- na.omit(gwas.markers)
gwas.markers <- stri_list2matrix(gwas.markers, byrow = FALSE)
gwas.markers <- as.data.frame(gwas.markers)
colnames(gwas.markers) <- gene.names

#GBJ
nrow(na.omit(gwas.pvalues$Sobic.001G000800))

sum(gwas.pvalues[,1])
while(colsum(gwas.pvalues) > 1){
  
}



