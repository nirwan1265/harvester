<<<<<<< HEAD
=======
## Run PCA
j <- 1
for(i in paste0("snpset.id",sprintf("%02d",1:10))){
  d = get(i)
  assign(paste0("pca",sprintf("%02d",j)), snpgdsPCA(get(paste0("gdsfile",sprintf("%02d",j))), snp.id = d, num.thread = 2))
  assign(i,d)
  j = j+1
}

>>>>>>> 25da95751769b8569fd0f5f8a24c74e2faa70901
#PCA analysis
data("hapmap_geno")

#Genotype
genotype <- as.data.frame(hapmap_geno$genotype)
head(genotype)
typeof(genotype) # list
class(genotype) # dataframe
#Tables with numeric genotype

sample.id <- hapmap_geno$sample.id
typeof(sample.id) # character
class(sample.id) #Character
sample.id
#Name of Samples or IDs
#Column name
<<<<<<< HEAD
#[1] "NA19152"     "NA19139"     "NA18912"     "NA19160"     "NA07034"    
#[6] "NA07055"     "NA12814"     "NA10847"     "NA18532"     "NA18561" 
=======

>>>>>>> 25da95751769b8569fd0f5f8a24c74e2faa70901

snp.id <- hapmap_geno$snp.id
typeof(snp.id) # character
class(snp.id) #Character
snp.id
#Row name
<<<<<<< HEAD
#[656] "rs12288829" "rs3133395"  "rs11222619" "rs7938283"  "rs1382840" 
#[661] "rs2239153"  "rs710415"   "rs10772627" "rs7311774"  "rs10846382"
=======
>>>>>>> 25da95751769b8569fd0f5f8a24c74e2faa70901

snp.position <- hapmap_geno$snp.position
typeof(snp.position) #double
class(snp.position) # integer
snp.position
#SNP position 
<<<<<<< HEAD
#[847]  11223058  11283630  12537074  13265827  14225059  14781966
#[853]  17468736  18121270  19639950  26516549  26652675  28767579
=======

>>>>>>> 25da95751769b8569fd0f5f8a24c74e2faa70901

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
####PCA


#Packages used:
library(vcfR)

#Genotype information table for running PCA
#Reading the path of the VCF files
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Data for sorghum/sorghum/Lasky.hapmap")
for(i in sprintf("%02d", 1:10)){
  assign(paste0("vcf.fn",i),paste0("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Data for sorghum/sorghum/Lasky.hapmap/hapmap.chr",i,".vcf"))
}

#Reading the VCF files:
j <- 1
for(i in paste0("vcf.fn", sprintf("%02d", 1:10))){
  d = get(i)
  assign(paste0("pca.geno.info",sprintf("%02d", j)), read.vcfR(d, verbose = FALSE))
  assign(i,d)
  j <- j + 1
}
#Save RDC files for VCF files. Takes long to load plus very large files
#j <- 1
#for(i in paste0("pca.geno.info", sprintf("%02d", 1:10))){
#  d = get(i)
#  saveRDS(d, paste0("pca.geno.info.chr",sprintf("%02d" , j),".RDS"))
#  assign(i,d)
#  j <- j+1
#}

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
for(i in sprintf("%02d", 1:10)){
  assign(paste0("geno", i), readRDS(paste0("geno",i,".RDS")))
}

#Collecting all the names
sample.id <- tidyr::gather(gwas01.gene.names)
sample.id <- sample.id[,-1]
sample.id <- as.data.frame(na.omit(sample.id))

#Collecting all the SNPs
library(tidyr)
x <- tidyr::gather(gwas01.Marker)
x <- x[,-1]
x <- as.data.frame(na.omit(x))
x <- as.vector(x$`na.omit(x)`)
x

#Sub-setting only the required SNPs
xy <- geno01[x]
xy[] <- lapply(xy, as.integer)

#Sub-setting the required SNPs from snp informations
x <- as.data.frame(x)
colnames(x) <- "ID"
xyz <- inner_join(pca.snp.info01,x, by = "ID")



#Convert to GDS
xyza <- snpgdsCreateGeno("test.gds",genmat = xy,
                 sample.id = xyz$ID,
                 snp.id = xyz$ID,
                 snp.chromosome = snp.chromosome01,
                 snp.position = snp.position01,
                 snp.allele = snp.allele01,
                 snpfirstdim = TRUE)





install.packages("vcfR")
library(vcfR)

chr1.vcf <- read.vcfR(vcf.fn01, verbose = FALSE)
x <- chr1.vcf@gt
















#SNP.id
j <- 1
for(i in paste0("pca.snp.info", sprintf("%02d", 1:10))){
  d = get(i)
  assign(paste0("snp.id",sprintf("%02d", j)), d$ID)
  assign(i,d)
  j <- j + 1
}

#SNP.position
j <- 1
for(i in paste0("pca.snp.info", sprintf("%02d", 1:10))){
  d = get(i)
  assign(paste0("snp.position",sprintf("%02d", j)), d$POS)
  assign(i,d)
  j <- j + 1
}
for(i in paste0("snp.position", sprintf("%02d", 1:10))){
  d = get(i)
  d <- as.double(d)
  assign(i,d)
}

#SNP.allele
j <- 1
for(i in paste0("pca.snp.info", sprintf("%02d", 1:10))){
  d = get(i)
  assign(paste0("snp.allele",sprintf("%02d", j)), paste(d$REF,d$ALT, sep = "/"))
  assign(i,d)
  j <- j + 1
}

#SNP.chromosome
j <- 1
for(i in paste0("pca.snp.info", sprintf("%02d", 1:10))){
  d = get(i)
  assign(paste0("snp.chromosome",sprintf("%02d", j)), d$CHROM)
  assign(i,d)
  j <- j + 1
}
for(i in paste0("snp.chromosome", sprintf("%02d", 1:10))){
  d = get(i)
  d <- as.integer(d)
  assign(i,d)
}

j <- 1
for(i in paste0("pca.geno", sprintf("%02d", 1:10))){
  d = get(i)
  d = d[,-1]
  assign(i,d)
  j <- j + 1
}

#Sample.ID
sample.id <- "REFERENCE_GENOME"
typeof(sample.id)
