###Packages required
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

###Local Functions:
##Gene Name filtering
split.names <- function(x,split){
  split.genename <- unlist(strsplit(x, split = ';', fixed = TRUE))[2]
  split.genename2 <- unlist(strsplit(split.genename, split = "=", fixed = TRUE))[2]
  return(split.genename2)
}

##Pvalues Combinations Test
pvalue.combine <- function(gwas.fstat, gwas.markers, gwas.pvalue, geno, tab.pc){
  x <- as.data.frame(matrix(0, nrow = 1, ncol = 1))
  y <- vector()
  combined.test.statistics <- as.data.frame(matrix(NA, nrow = 1, ncol = 1))
  
  for (i in 1:ncol(gwas.fstat)){
    for(j in 1:sum(!is.na(gwas.fstat[,i]))){
      x[1,j] <- as.double(gwas.fstat[j,i])
      y[j] <- as.vector(as.character(gwas.markers[j,i]))
    }
    if(ncol(x) >= 2){
      x <- as.matrix(as.double(x))
      ref_genotype <- as.data.frame(geno10[,colnames(geno10) %in% y])
      cor_mat <- estimate_ss_cor(ref_pcs=tab.pc, ref_genotypes=ref_genotype, link_function='log')
      bj.test <- BJ(test_stats = x, cor_mat=cor_mat)
      gbj.test <- GBJ(test_stats = x, cor_mat=cor_mat)
      minP.test <- minP(test_stats = x, cor_mat=cor_mat)
      hc.test <- HC(test_stats = x, cor_mat=cor_mat)
      ghc.test <- GHC(test_stats = x, cor_mat=cor_mat)
      combined.test.statistics[i,1] <- bj.test$BJ_pvalue
      combined.test.statistics[i,2] <- gbj.test$GBJ_pvalue
      combined.test.statistics[i,3] <- hc.test$HC_pvalue
      combined.test.statistics[i,4] <- ghc.test$GHC_pvalue
      combined.test.statistics[i,5] <- minP.test$minP_pvalue
      x <- as.data.frame(matrix(0, nrow = 1, ncol = 1))
      y <- vector()
    }else{
      combined.test.statistics[i,1] <- as.double(gwas.pvalue[1,i])
      combined.test.statistics[i,2] <- as.double(gwas.pvalue[1,i])
      combined.test.statistics[i,3] <- as.double(gwas.pvalue[1,i])
      combined.test.statistics[i,4] <- as.double(gwas.pvalue[1,i])
      combined.test.statistics[i,5] <- as.double(gwas.pvalue[1,i])
    }
  }
  
}


###Loading database
##Download SNP database
#SNP Database 
#http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-52/plants/gff3/sorghum_bicolor/
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Data for sorghum/sorghum/SNP annotation/Sorghum")
snp.db <- read.table(file ="Sbicolor_454_v3.1.1.gene_exons.gff3", sep = "\t", header = FALSE)
colnames(snp.db) <- c("Chromosome","Database","Region","Start","End","NA","Strand","NA2","Gene")

###Loading query SNPs for GWAS from RDS file
##Location in the server: /rsstu/users/r/rrellan/sara/SorghumGEA/results/GLM_20220222 - GLM
##Location in the server: /rsstu/users/r/rrellan/sara/SorghumGEA/results/GLM_20220224 - GLM
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Data for sorghum/sorghum/GWAS.results")
for(i in sprintf("%02d", 1:10)){
  assign(paste0("query.snp.gwas", i) , readRDS(file = paste0("glm_sol_VL_",i,"_20220224_03_09.RDS")))
}

for(i in sprintf("%02d", 1:10)){
  assign(paste0("query.snp.gwas",i), data.frame(get(paste0("query.snp.gwas",i))[1]))
  assign(paste0("query.snp.gwas",i), get(paste0("query.snp.gwas",i))[,c(2,3,4,5,6)])
}


###Sub-setting chromosomes from db
for(i in sprintf("%02d", 1:10)){
  assign(paste0("x",i), paste0("Chr",i))
}
for(i in sprintf("%02d", 1:10)){
  assign(paste0("db.",i), snp.db[which(snp.db$Chromosome == get(paste0("x",i))), ])
}


##Making GRanges for Database 
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Data for sorghum/sorghum/GWAS.results")
for(i in sprintf("%02d", 1:10)){
  assign(paste0("gr.db", i) , GRanges(seqnames = paste0("chr",i), ranges = IRanges(start = get(paste0("db.",i))[,"Start"], end = get(paste0("db.",i))[,"End"]), strand = get(paste0("db.",i))[,"Strand"], Region = get(paste0("db.",i))[,"Region"], Gene = get(paste0("db.",i))[,"Gene"]))
}

##Making GRanges for gwas Query
for(i in sprintf("%02d", 1:10)){
  assign(paste0("gr.q", i) , GRanges(seqnames = paste0("chr",i), ranges = IRanges(start = get(paste0("query.snp.gwas",i))[,"GLM_Stats.Pos"], width =1, fstat = get(paste0("query.snp.gwas",i))[,"GLM_Stats.marker_F"], Marker = get(paste0("query.snp.gwas",i))[,"GLM_Stats.Marker"],pvalue = get(paste0("query.snp.gwas",i))[,"GLM_Stats.p"])))
}

##Overlaps
for(i in sprintf("%02d", 1:10)){
  assign(paste0("common",i), as.data.frame(findOverlapPairs(get(paste0("gr.db",i)), get(paste0("gr.q",i)))))
}


###Combining F-stat:
##Filtering out the gene, fstat, Marker and pvalue columns
#Making a new dummy dataframe
for(i in sprintf("%02d", 1:10)){
  assign(paste0("gwas",i), get(paste0("common",i)))
}

#Subsetting only gene
gwas01 <- gwas01[which(gwas01$first.X.Region == "gene"), ]
gwas02 <- gwas02[which(gwas02$first.X.Region == "gene"), ]
gwas03 <- gwas03[which(gwas03$first.X.Region == "gene"), ]
gwas04 <- gwas04[which(gwas04$first.X.Region == "gene"), ]
gwas05 <- gwas05[which(gwas05$first.X.Region == "gene"), ]
gwas06 <- gwas06[which(gwas06$first.X.Region == "gene"), ]
gwas07 <- gwas07[which(gwas07$first.X.Region == "gene"), ]
gwas08 <- gwas08[which(gwas08$first.X.Region == "gene"), ]
gwas09 <- gwas09[which(gwas09$first.X.Region == "gene"), ]
gwas10 <- gwas10[which(gwas10$first.X.Region == "gene"), ]
#Filtering out the required columns
for(i in sprintf("%02d", 1:10)){
  assign(paste0("gwas",i), data.frame(get(paste0("gwas",i)))[,c(7,15,16,17)])
}
#Renaming columns
for(i in paste0("gwas", sprintf("%02d", 1:10))){
  d=get(i)
  colnames(d) = c("Gene","fstat","Marker","pvalue")
  assign(i,d)
}

##Sorting a/c gene name
#Mayber not required, but might be for other database or when combining all result
for(i in paste0("gwas", sprintf("%02d", 1:10))){
  d=get(i)
  d <- d[gtools::mixedorder(d$Gene), ]
  assign(i,d)
}

##Table with F-stat values
for(i in paste0("gwas",sprintf("%02d",1:10))){
  d=get(i)
  assign(paste0(i,".fstat"), dcast(setDT(d), Gene~rowid(Gene, prefix = "fstat"), value.var = "fstat"))
  assign(i,d)
}
#----> continue
##Adding gene names
#Table with all the gene names from all the chromosome
#Can use this table for all the other tables with Marker and pvalue data
for(i in paste0("gwas", sprintf("%02d", 1:10))){
  d=get(i)
  assign(paste0(i,".gene.names"), get(paste0(i,".fstat"))[,1])
  assign(i,d)
}
for(i in paste0("gwas", 10,".gene.names")){
  d=get(i)
  d <-  as.data.frame(t(apply(d,1,split.names)))
  assign(i,d)
}

#----> continue
for(i in 1:10){
  assign(paste0("gwas.fstat",i), dcast(setDT(paste0("gwas",i)), Gene~rowid(Gene, prefix = "fstat"), value.var = "fstat"))
}
#Adding gene names
for(i in 1:10){
  assign(paste0("gene.names",i), paste0("gwas.fstat",i)[1, ])
}
for(i in 1:10){
  assign(paste0("gene.names",i), apply(paste0("gene.names",i), 2, split.names))
  assign(paste0("gene.names",i), as.data.frame(t(paste0("gene.names"),i)))
}
for(i in 1:10){
  assign(paste0("gwas.fstat",i), paste0("gwasfstat",i)[-1, ])
  assign(paste0("gwas.fstat",i), apply(paste0("gwasfstat",i),2,sort))
  assign(paste0("gwas.fstat",i), na.omit(paste0("gwasfstat",i)))
  assign(paste0("gwas.fstat",i), stri_list2matrix(paste0("gwasfstat",i), byrow = FALSE))
  colnames(assign(paste0("gwas.fstat",i), paste0("gene.names",i)))
}



##Table with SNP markers
for(i in 10){
  assign(paste0("gwas.markers",i), dcast(setDT(paste0("gwas",i)), Gene~rowid(Gene, prefix = "markers"), value.var = "markers"))
}

for(i in 1:10){
  assign(paste0("gwas.markers",i), dcast(setDT(paste0("gwas",i)), Gene~rowid(Gene, prefix = "markers"), value.var = "markers"))
}
for(i in 1:10){
  assign(paste0("gwas.markers",i), paste0("gwasmarkers",i)[-1, ])
  assign(paste0("gwas.markers",i), apply(paste0("gwasmarkers",i),2,sort))
  assign(paste0("gwas.markers",i), na.omit(paste0("gwasmarkers",i)))
  assign(paste0("gwas.markers",i), stri_list2matrix(paste0("gwasmarkers",i), byrow = FALSE))
  colnames(assign(paste0("gwas.markers",i), paste0("gene.names",i)))
}


##Table with pvalues

for(i in 10){
  assign(paste0("gwas.pvalue",i), dcast(setDT(paste0("gwas",i)), Gene~rowid(Gene, prefix = "pvalue"), value.var = "pvalue"))
}

for(i in 1:10){
  assign(paste0("gwas.pvalue",i), dcast(setDT(paste0("gwas",i)), Gene~rowid(Gene, prefix = "pvalue"), value.var = "pvalue"))
}
for(i in 1:10){
  assign(paste0("gwas.pvalue",i), paste0("gwaspvalue",i)[-1, ])
  assign(paste0("gwas.pvalue",i), apply(paste0("gwaspvalue",i),2,sort))
  assign(paste0("gwas.pvalue",i), na.omit(paste0("gwaspvalue",i)))
  assign(paste0("gwas.pvalue",i), stri_list2matrix(paste0("gwaspvalue",i), byrow = FALSE))
  colnames(assign(paste0("gwas.pvalue",i), paste0("gene.names",i)))
}


###PCA analysis required for GBJ
#Need a vcf file format of the hapmap which is converted to GDS format
#TASSEL or PLINK is used for converting hapmap to VCF file format
#Need a directory to  create the gds file. If working on the server, we might need to define this before starting

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Data for sorghum/sorghum/Sorghum.genotype")
for(i in 10){
  assign(paste0("vcf.fn",i),"path/to/file/","sorghum.chr",i,".vcf")
}

for(i in 1:10){
  assign(paste0("vcf.fn",i),"path/to/file/","sorghum.chr",i,".vcf")
}

For(i in 1:10){
  snpgdsVCF2GDS(paste0("vcf.fn",i), paste0("sorghum.ch",i,".gds"), method = "biallelic.only" , )
}

##Get GDS file data
for(i in 1:10){
  assign(paste0("genofile",i), snpgdsOpen(paste0("sorghum.ch",i,".gds")))
}

##LD-based SNP pruning
set.seed(1000)
# Try different LD thresholds for sensitivity analysis but read in a paper somewhere that 0.2 was used for GBJ
for(i in 10){
  assign(paste0("snpset",i), snpgdsLDpruning(paste0("genofile",i), ld.threshold = 0.2))
}
## Get all selected snp id
for(i in 1:10){
  assign(paste0("snpset",i), snpgdsLDpruning(paste0("genofile",i), ld.threshold = 0.2))
}
## Get all selected snp id
for(i in 1:10){
  assign(paste0("snpset.id",i), unlist(unname(paste0("snpset",i))))
}

## Run PCA
for(i in 1:10){
  assign(paste0("pca",i), snpgdsPCA(paste0("genofile",i), snp.id = paste0("snpset.id",i), num.thread = 2))
}


#In case there are population information
#https://www.bioconductor.org/packages/devel/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html
#In the case of no prior population information,
#Make a data.frame of eigen values
for(i in 1:10){
  assign(paste0("tab",i), data.frame(paste0("sample.id",i) = paste0("pca$sample.id",i),
                                     EV1 = paste0("pca",i,"$eigenvect")[,1],
                                     EV2 = paste0("pca",i,"$eigenvect")[,2],
                                     EV3 = paste0("pca",i,"$eigenvect")[,3],
                                     EV4 = paste0("pca",i,"$eigenvect")[,4],
                                     EV5 = paste0("pca",i,"$eigenvect")[,5],
                                     stringAsFactors = FALSE))
}

###Pvalue combination
##Pre-processing for pvalue combination
for(i in 1:10){
  assign(paste0("tab.pc",1), paste0("tab",i)[,c(2:6)])
}

##Numerical hapmap genotype file
#First need to impute the hapmap file and then change to numeric format and vice-versa
#IMPORTANT: This can be done in TASSEL. Remove <marker> and <numerical> text from the transposed txt file before loading 
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Data for sorghum/sorghum/Sorghum.genotype")
for(i in 1:10){
  assign(paste0("geno",i), read.table(file = paste0("numerical.imputed.hapmap",i,".txt"), header = TRUE, sep = "\t"))
}

##Combination tests using GBJ
for(i in 1:10){
  assign(paste0("pvalue.combine",i), pvalue.combine(paste0("gwas.fstat",i), paste0("gwas.markers",i), paste0("gwas.pvalue",i), paste0("geno",i),paste0("tab.pc",i)))
}

##Saving the result as RDS
for(i in 1:10){
  saveRDS(paste0("pvalue.combine",i), file = paste0("pvalue.combine.sorghum.chr",i))
}


