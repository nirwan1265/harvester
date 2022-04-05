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
#Tables with numeric genotype

sample.id <- hapmap_geno$sample.id
typeof(sample.id) # character
class(sample.id) #Character
sample.id
#Name of Samples or IDs
#Column name


snp.id <- hapmap_geno$snp.id
typeof(snp.id) # character
class(snp.id) #Character
snp.id
#Row name

snp.position <- hapmap_geno$snp.position
typeof(snp.position) #double
class(snp.position) # integer
snp.position
#SNP position 


snp.allele <- hapmap_geno$snp.allele
typeof(snp.allele) #character
class(snp.allele) # character
snp.allele
#[1] "A/G" "C/T" "A/C" "A/G" "A/G" "C/T" "C/T" "C/T" "G/T" "A/G" "A/T" "A/G" "A/G" "A/G" "G/T" "C/T" "A/G" "A/G" "C/T"
#[20] "A/G" "C/T" "C/T" "A/T" "A/G" "C/T" "A/G" "C/T" "G/T" "C/T" "A/G" "C/T" "C/T" "C/T" "A/G" "A/G" "A/C" "A/G" "C/T"
#Both alleles


install.packages("vcfR")
library(vcfR)

chr1.vcf <- read.vcfR(vcf.fn01, verbose = FALSE)
x <- chr1.vcf@gt
