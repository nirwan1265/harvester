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


<<<<<<< HEAD
=======


>>>>>>> 5df7a0819959d81625ea8743ef11f389cb69e873
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


###Subsetting chr from db
for(i in sprintf("%02d", 1:10)){
  assign(paste0("db.",i), get(snp.db[which(snp.db$Chromosome == paste0('"Chr',i,'"'))]))
}
for(i in 10){
  assign(paste0("db.",i), get(cat(sprintf("snp.db[which(snp.db$Chromosome == \"Chr%02d\"), ]",i))))
}

##Making GRanges for Database 
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Data for sorghum/sorghum/GWAS.results")
for(i in 010){
  assign(paste0("gr.db", i) , GRanges(seqnames = paste0("chr",i), ranges = IRanges(start = get(paste0("db.",i))[,"Start"], end = get(paste0("db.",i))[,"End"]), strand = get(paste0("db.",i))[,"Strand"], Region = get(paste0("db.",i))[,"Region"], Gene = get(paste0("db.",i))[,"Gene"]))
}

#Making GRanges for gwas Query
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Data for sorghum/sorghum/GWAS.results")
for(i in 10){
  assign(paste0("query.snp.gwas",i), readRDS(file = paste0("glm_sol_VL_",i,"_20220222_13_55.RDS")))
}
for(i in 10){
  assign(paste0("query.snp.gwas",i), as.data.frame(paste0("query.snp.gwas",i,"$GLM_Stats")))
  assign(paste0("query.snp.gwas",i), paste0("query.snp.gwas",i[,c(2,3,4,5,6)]))
}

for(i in 10){
  assign(paste0("gr.q", i) , GRanges(seqnames = paste0("chr",i), ranges = IRanges(start = get(paste0("query",i,".gwas"))[,"Pos"], width =1, fstat = get(paste0("query",i,".gwas"))[,"p"], maker = get(paste0("query",i,".gwas"))[,"Marker"])))
}


##Overlaps
for(i in 10){
  assign(paste0("common",i), as.data.frame(findOverlapPairs(paste0("gr.db",i), paste0("gr.q.",i))))
}


###Combining F-stat:
##Filtering out only the gene and fstat columns
<<<<<<< HEAD
for(i in 10){
  assign(paste0("gwas",i), paste0("common",i))
}
for(i in 10){
=======
for(i in 1:10){
  assign(paste0("gwas",i), paste0("common",i))
}
for(i in 1:10){
>>>>>>> 5df7a0819959d81625ea8743ef11f389cb69e873
  assign(paste0("gwas",i), paste0("gwas",i)[,c(6,7,15,16,17)])
  assign(paste0("gwas",i), paste0("gwas",i)[which(paste0("gwas",i,"first.X.Region") == "gene"), ])
  assign(paste0("gwas",i), paste0("gwas",i)[,c(2,3,4,5)])
}
<<<<<<< HEAD
for(i in 10){
=======
for(i in 1:10){
>>>>>>> 5df7a0819959d81625ea8743ef11f389cb69e873
  colnames(assign(paste0("gwas",i)), c("Gene","fstat","Marker","pvalue"))
}

##Sorting a/c gene name
<<<<<<< HEAD
for(i in 10){
=======
for(i in 1:10){
>>>>>>> 5df7a0819959d81625ea8743ef11f389cb69e873
  assign(paste0("gwas",i), paste0("gwas",i)[gtools::mixedorder(paste0("gwas",i,"Gene")), ])
}

##Table with F-stat values
<<<<<<< HEAD
for(i in 10){
  assign(paste0("gwas.fstat",i), dcast(setDT(paste0("gwas",i)), Gene~rowid(Gene, prefix = "fstat"), value.var = "fstat"))
}
#Adding gene names
for(i in 10){
  assign(paste0("gene.names",i), paste0("gwas.fstat",i)[1, ])
}
for(i in 10){
  assign(paste0("gene.names",i), apply(paste0("gene.names",i), 2, split.names))
  assign(paste0("gene.names",i), as.data.frame(t(paste0("gene.names"),i)))
}
for(i in 10){
=======
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
>>>>>>> 5df7a0819959d81625ea8743ef11f389cb69e873
  assign(paste0("gwas.fstat",i), paste0("gwasfstat",i)[-1, ])
  assign(paste0("gwas.fstat",i), apply(paste0("gwasfstat",i),2,sort))
  assign(paste0("gwas.fstat",i), na.omit(paste0("gwasfstat",i)))
  assign(paste0("gwas.fstat",i), stri_list2matrix(paste0("gwasfstat",i), byrow = FALSE))
  colnames(assign(paste0("gwas.fstat",i), paste0("gene.names",i)))
}



##Table with SNP markers
<<<<<<< HEAD
for(i in 10){
  assign(paste0("gwas.markers",i), dcast(setDT(paste0("gwas",i)), Gene~rowid(Gene, prefix = "markers"), value.var = "markers"))
}
for(i in 10){
=======
for(i in 1:10){
  assign(paste0("gwas.markers",i), dcast(setDT(paste0("gwas",i)), Gene~rowid(Gene, prefix = "markers"), value.var = "markers"))
}
for(i in 1:10){
>>>>>>> 5df7a0819959d81625ea8743ef11f389cb69e873
  assign(paste0("gwas.markers",i), paste0("gwasmarkers",i)[-1, ])
  assign(paste0("gwas.markers",i), apply(paste0("gwasmarkers",i),2,sort))
  assign(paste0("gwas.markers",i), na.omit(paste0("gwasmarkers",i)))
  assign(paste0("gwas.markers",i), stri_list2matrix(paste0("gwasmarkers",i), byrow = FALSE))
  colnames(assign(paste0("gwas.markers",i), paste0("gene.names",i)))
}


##Table with pvalues
<<<<<<< HEAD
for(i in 10){
  assign(paste0("gwas.pvalue",i), dcast(setDT(paste0("gwas",i)), Gene~rowid(Gene, prefix = "pvalue"), value.var = "pvalue"))
}
for(i in 10){
=======
for(i in 1:10){
  assign(paste0("gwas.pvalue",i), dcast(setDT(paste0("gwas",i)), Gene~rowid(Gene, prefix = "pvalue"), value.var = "pvalue"))
}
for(i in 1:10){
>>>>>>> 5df7a0819959d81625ea8743ef11f389cb69e873
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
<<<<<<< HEAD
for(i in 10){
  assign(paste0("vcf.fn",i),"path/to/file/","sorghum.chr",i,".vcf")
}

For(i in 10){
=======
for(i in 1:10){
  assign(paste0("vcf.fn",i),"path/to/file/","sorghum.chr",i,".vcf")
}

For(i in 1:10){
>>>>>>> 5df7a0819959d81625ea8743ef11f389cb69e873
  snpgdsVCF2GDS(paste0("vcf.fn",i), paste0("sorghum.ch",i,".gds"), method = "biallelic.only" , )
}

##Get GDS file data
<<<<<<< HEAD
for(i in 10){
=======
for(i in 1:10){
>>>>>>> 5df7a0819959d81625ea8743ef11f389cb69e873
  assign(paste0("genofile",i), snpgdsOpen(paste0("sorghum.ch",i,".gds")))
}

##LD-based SNP pruning
set.seed(1000)
# Try different LD thresholds for sensitivity analysis but read in a paper somewhere that 0.2 was used for GBJ
<<<<<<< HEAD
for(i in 10){
  assign(paste0("snpset",i), snpgdsLDpruning(paste0("genofile",i), ld.threshold = 0.2))
}
## Get all selected snp id
for(i in 10){
=======
for(i in 1:10){
  assign(paste0("snpset",i), snpgdsLDpruning(paste0("genofile",i), ld.threshold = 0.2))
}
## Get all selected snp id
for(i in 1:10){
>>>>>>> 5df7a0819959d81625ea8743ef11f389cb69e873
  assign(paste0("snpset.id",i), unlist(unname(paste0("snpset",i))))
}

## Run PCA
<<<<<<< HEAD
for(i in 10){
=======
for(i in 1:10){
>>>>>>> 5df7a0819959d81625ea8743ef11f389cb69e873
  assign(paste0("pca",i), snpgdsPCA(paste0("genofile",i), snp.id = paste0("snpset.id",i), num.thread = 2))
}


#In case there are population information
#https://www.bioconductor.org/packages/devel/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html
#In the case of no prior population information,
#Make a data.frame of eigen values
<<<<<<< HEAD
for(i in 10){
=======
for(i in 1:10){
>>>>>>> 5df7a0819959d81625ea8743ef11f389cb69e873
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
<<<<<<< HEAD
for(i in 10){
=======
for(i in 1:10){
>>>>>>> 5df7a0819959d81625ea8743ef11f389cb69e873
  assign(paste0("tab.pc",1), paste0("tab",i)[,c(2:6)])
}

##Numerical hapmap genotype file
#First need to impute the hapmap file and then change to numeric format and vice-versa
#IMPORTANT: This can be done in TASSEL. Remove <marker> and <numerical> text from the transposed txt file before loading 
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Data for sorghum/sorghum/Sorghum.genotype")
<<<<<<< HEAD
for(i in 10){
=======
for(i in 1:10){
>>>>>>> 5df7a0819959d81625ea8743ef11f389cb69e873
  assign(paste0("geno",i), read.table(file = paste0("numerical.imputed.hapmap",i,".txt"), header = TRUE, sep = "\t"))
}

##Combination tests using GBJ
<<<<<<< HEAD
for(i in 10){
=======
for(i in 1:10){
>>>>>>> 5df7a0819959d81625ea8743ef11f389cb69e873
  assign(paste0("pvalue.combine",i), pvalue.combine(paste0("gwas.fstat",i), paste0("gwas.markers",i), paste0("gwas.pvalue",i), paste0("geno",i),paste0("tab.pc",i)))
}

##Saving the result as RDS
<<<<<<< HEAD
for(i in 10){
=======
for(i in 1:10){
>>>>>>> 5df7a0819959d81625ea8743ef11f389cb69e873
  saveRDS(paste0("pvalue.combine",i), file = paste0("pvalue.combine.sorghum.chr",i))
}


