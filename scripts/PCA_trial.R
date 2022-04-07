## Run PCA
j <- 1
for(i in paste0("snpset.id",sprintf("%02d",1:10))){
  d = get(i)
  assign(paste0("pca",sprintf("%02d",j)), snpgdsPCA(get(paste0("gdsfile",sprintf("%02d",j))), snp.id = d, num.thread = 2))
  assign(i,d)
  j = j+1
}

#PCA analysis
data("hapmap_geno")

#Genotype
genotype <- as.data.frame(hapmap_geno$genotype)
head(genotype)
typeof(genotype) # list
class(genotype) # dataframe
genotype



#Tables with numeric genotype
sample.id <- hapmap_geno$sample.id
typeof(sample.id) # character
class(sample.id) #Character
sample.id

#Name of Samples or IDs
#Column name
#[1] "NA19152"     "NA19139"     "NA18912"     "NA19160"     "NA07034"    
#[6] "NA07055"     "NA12814"     "NA10847"     "NA18532"     "NA18561" 

snp.id <- hapmap_geno$snp.id
typeof(snp.id) # character
class(snp.id) #Character
snp.id
#Row name
#[656] "rs12288829" "rs3133395"  "rs11222619" "rs7938283"  "rs1382840" 
#[661] "rs2239153"  "rs710415"   "rs10772627" "rs7311774"  "rs10846382"

snp.position <- hapmap_geno$snp.position
typeof(snp.position) #double
class(snp.position) # integer
snp.position
#SNP position 
#[847]  11223058  11283630  12537074  13265827  14225059  14781966
#[853]  17468736  18121270  19639950  26516549  26652675  28767579


snp.allele <- hapmap_geno$snp.allele
typeof(snp.allele) #character
class(snp.allele) # character
snp.allele
#[1] "A/G" "C/T" "A/C" "A/G" "A/G" "C/T" "C/T" "C/T" "G/T" "A/G" "A/T" "A/G" "A/G" "A/G" "G/T" "C/T" "A/G" "A/G" "C/T"
#[20] "A/G" "C/T" "C/T" "A/T" "A/G" "C/T" "A/G" "C/T" "G/T" "C/T" "A/G" "C/T" "C/T" "C/T" "A/G" "A/G" "A/C" "A/G" "C/T"
#Both alleles


snp.chromosome <- hapmap_geno$snp.chromosome
typeof(snp.chromosome) #integer
class(snp.chromosome) #integer
snp.chromosome
# [1]  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1


####PCA########################################################################
#Genotype information table for running PCA
#Reading the path of the VCF files
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Data for sorghum/sorghum/Lasky.hapmap")
for(i in sprintf("%02d", 1:10)){
  assign(paste0("vcf.fn",i),paste0("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Data for sorghum/sorghum/Lasky.hapmap/hapmap.chr",i,".vcf"))
}

#Reading the VCF files:
#j <- 1
#for(i in paste0("vcf.fn", sprintf("%02d", 1:10))){
#  d = get(i)
#  assign(paste0("pca.geno.info",sprintf("%02d", j)), read.vcfR(d, verbose = FALSE))
#  assign(i,d)
#  j <- j + 1
#}
#Save RDC files for VCF files. Takes long to load plus very large files
#j <- 1
#for(i in paste0("pca.geno.info", sprintf("%02d", 1:10))){
#  d = get(i)
#  saveRDS(d, paste0("pca.geno.info.chr",sprintf("%02d" , j),".RDS"))
#  assign(i,d)
#  j <- j+1
#}

#Reading PCA genotype RDS file
for(i in sprintf("%02d", 1:10)){
  assign(paste0("pca.geno.info", i) , readRDS(file = paste0("pca.geno.info.chr",i,".RDS")))
}


#Getting SNPs info
j <- 1
for(i in paste0("pca.geno.info", sprintf("%02d", 1:10))){
  d = get(i)
  assign(paste0("pca.snp.info",sprintf("%02d", j)), as.data.frame(d@fix))
  assign(i,d)
  j <- j + 1
}
for(i in paste0("pca.snp.info", sprintf("%02d", 1:10))){
  d = get(i)
  d <- d[,c(1:5)]
  assign(i,d)
}
for(i in paste0("pca.snp.info", sprintf("%02d", 1:10))){
  d <-  get(i)
  d$allele <- paste0(d$REF,"/", d$ALT) 
  assign(i,d)
}

#Read numerical genotype 
for(i in sprintf("%02d", 1)){
  assign(paste0("geno", i), readRDS(paste0("geno",i,".RDS")))
}
geno01[1:50,1:50]

#Copying the genotype file in dummy variable
genofile <- as.data.frame(geno01)
genofile[1:6,1:6]

#Sub-setting only the required SNPs
sample.id <- genofile[,1]

#Removing the 

#Collecting all the SNPs
snp.collect <- tidyr::gather(gwas01.Marker)
snp.collect <- snp.collect[,-1]
snp.collect <- as.data.frame(na.omit(snp.collect))
snp.collect <- as.vector(snp.collect$`na.omit(snp.collect)`)
snp.collect

#Sub-setting the required SNPs from  SNP information
snp.collect <- as.data.frame(snp.collect)
colnames(snp.collect) <- "ID"
pca.snp.gds <- inner_join(pca.snp.info01,snp.collect, by = "ID")

#Converting Pos and CHR to numeric
pca.snp.gds$position <- as.integer(pca.snp.gds$POS)
pca.snp.gds$chromosome <- as.integer(pca.snp.gds$CHROM)

#Sub-setting the required SNPs from genotype file
snp.collect.t <- as.vector(t(snp.collect))
genofile.unique <- genofile[snp.collect.t]
genofile.unique <- as.matrix(t(genofile.unique))

#Convert to GDS
gdsfile <- snpgdsCreateGeno("test.gds",genmat = genofile.unique,
                 sample.id = sample.id,
                 snp.id = pca.snp.gds$ID,
                 snp.chromosome = pca.snp.gds$chromosome,
                 snp.position = pca.snp.gds$position,
                 snp.allele = pca.snp.gds$allele,
                 snpfirstdim = TRUE)


duplicated(pca.snp.gds$ID)



library(AssocTests)
data("drS.eg")

eigenstratG.eg <- matrix(rbinom(3000, 2, 0.5), ncol = 30)
write.table(eigenstratG.eg, file = "eigenstratG.eg.txt", quote = FALSE,
            sep = "", row.names = FALSE, col.names = FALSE)
x <- eigenstrat(genoFile = "eigenstratG.eg.txt", outFile.Robj = "eigenstrat.result.list",
                outFile.txt = "eigenstrat.result.txt", rm.marker.index = NULL,
                rm.subject.index = NULL, miss.val = 9, num.splits = 10,
                topK = NULL, signt.eigen.level = 0.01, signal.outlier = FALSE,
                iter.outlier = 5, sigma.thresh = 6)
file.remove("eigenstratG.eg.txt", "eigenstrat.result.list", "eigenstrat.result.txt")



#The genotype file
genofile <- as.data.frame(geno01)
genofile[1:6,1:6]


#Collecting all the Markers
library(tidyr)
x <- tidyr::gather(gwas01.Marker)
x <- x[,-1]
x <- as.data.frame(na.omit(x))
x <- as.vector(x$`na.omit(x)`)
head(x,6)

#Sub-setting only the required SNPs
genofile2 <- genofile[x]
genofile2 <- genofile2[-1,1:40]
#Replacing 0.5 with 9
genofile2[genofile2 == 0.5] <- 9


write.table(genofile2, file = "eigenstratG.eg.txt", quote = FALSE,
            sep = "", row.names = FALSE, col.names = FALSE)
x <- eigenstrat(genoFile = "eigenstratG.eg.txt", outFile.Robj = "eigenstrat.result.list",
                outFile.txt = "eigenstrat.result.txt", rm.marker.index = NULL,
                rm.subject.index = NULL, miss.val = 9, num.splits = 10,
                topK = NULL, signt.eigen.level = 0.01, signal.outlier = FALSE,
                iter.outlier = 5, sigma.thresh = 6)
