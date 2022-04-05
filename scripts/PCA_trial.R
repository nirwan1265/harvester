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
####PCA


#Packages used:
library(vcfR)



#Reading the vcf file
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

#Getting SNPs info
j <- 1
for(i in paste0("pca.geno.info", sprintf("%02d", 1:10))){
  d = get(i)
  assign(paste0("pca.snp.info",sprintf("%02d", j)), as.data.frame(d@fix))
  assign(i,d)
  j <- j + 1
}


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

#Read numerical genotype 
for(i in sprintf("%02d", 1:10)){
  assign(paste0("pca.geno", i), read.table(paste0("PCA.chr",i,".txt"), header = TRUE))
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

#Convert to GDS
snpgdsCreateGeno("test.gds",genmat = pca.geno01,
                 sample.id = sample.id,
                 snp.id = snp.id01,
                 snp.chromosome = snp.chromosome01,
                 snp.position = snp.position01,
                 snp.allele = snp.allele01,
                 snpfirstdim = TRUE)





