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
#for(i in sprintf("%02d", 1:10)){
#  assign(paste0("vcf.fn",i),paste0("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Data for sorghum/sorghum/Lasky.hapmap/hapmap.chr",i,".vcf"))
#}

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
for(i in sprintf("%02d", 10)){
  assign(paste0("geno", i), readRDS(paste0("geno",i,".RDS")))
}
geno10[1:6,1:6]


#Copying the genotype file in dummy variable
genofile <- as.data.frame(geno10)
genofile[1:6,1:6]

#Sample.id
sample.id <- as.data.frame(genofile[,1])

#Collecting all the SNPs
snp.collect <- tidyr::gather(gwas01.Marker)
snp.collect <- snp.collect[,-1]
snp.collect <- as.data.frame(na.omit(snp.collect))
snp.collect <- as.vector(snp.collect$`na.omit(snp.collect)`)
snp.collect

#Sub-setting the required SNPs from  SNP information
snp.collect <- as.data.frame(snp.collect)
colnames(snp.collect) <- "ID"
pca.snp.gds <- right_join(pca.snp.info01,snp.collect, by = "ID")
#Removing duplicate values
pca.snp.gds <- pca.snp.gds[!duplicated(pca.snp.gds$position), ]

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



##EIGENSTRAT

data("drS.eg")

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Data for sorghum/sorghum/Lasky.hapmap")


library(AssocTests)
data("drS.eg")


#Read numerical genotype 
for(i in sprintf("%02d", 1)){
  assign(paste0("geno", i), readRDS(paste0("geno",i,".RDS")))
}
geno01[1:6,1:6]

#Collecting all the SNPs
snp.collect <- tidyr::gather(gwas01.Marker)
snp.collect <- snp.collect[,-1]
snp.collect <- as.data.frame(na.omit(snp.collect))
snp.collect <- as.vector(snp.collect$`na.omit(snp.collect)`)
snp.collect
snp.collect


#Subsetting the required SNPs
geno01.subset <- geno01[snp.collect]

#Replacing 0.5 with 9
geno01.subset[geno01.subset == 0.5] <- 9


#Subsetting only SNPs for genes
head(gwas01.Marker)
gwas01.Marker[1:6,1:6]
sum(!is.na(gwas01.Marker$Sobic.001G000100))

snp.tbl <- as.data.frame(matrix(NA))
geno.tbl <- as.data.frame(matrix(NA))
snp.tbl <- as.vector(gwas10.Marker[,2])
snp.tbl <- snp.tbl[!is.na(snp.tbl)]
geno.tbl <- geno10[snp.tbl]

geno.tbl[geno.tbl == 0.5] <- 9
=======
geno.tbl[geno.tbl == 9] <- 1
>>>>>>> add704465cad8d5bc3d2e2a38b769f5a90d80dac

system("pwd")

write.table(geno.tbl, file = "eigenstratG.eg.txt", quote = FALSE,
            sep = "", row.names = FALSE, col.names = FALSE)
x <- eigenstrat(genoFile = "eigenstratG.eg.txt", outFile.Robj = "eigenstrat.result.list",
                outFile.txt = "eigenstrat.result.txt", rm.marker.index = NULL,

                rm.subject.index = NULL, miss.val = 9, num.splits = 1,
                topK = NULL, signt.eigen.level = 0.01, signal.outlier = FALSE,
                iter.outlier = 5, sigma.thresh = 6)
x$eigenvectors

??eigenstrat

=======
                rm.subject.index = NULL, miss.val = 9, num.splits = 10,
                topK = NULL, signt.eigen.level = 0.01, signal.outlier = FALSE,
                iter.outlier = 5, sigma.thresh = 6)
>>>>>>> add704465cad8d5bc3d2e2a38b769f5a90d80dac
x <- x$eigenvectors
x <- as.matrix(x[,1:2])
x
typeof(x)
class(x)

FGFR2_cor_mat <- estimate_ss_cor(ref_pcs=x, ref_genotypes=geno.tbl, link_function='logit')


for(i in ncol(gwas01.Marker)){
  while(sum(length(which(!is.na(gwas01.Marker[,i])))) > 3){
    snp.tbl <- as.vector(gwas01.Marker[,i])
    snp.tbl <- snp.tbl[!is.na(snp.tbl)]
    geno.tbl <- geno01[snp.tbl]
    write.table(geno.tbl, file = "eigenstratG.eg.txt", quote = FALSE,

                sep = " ", row.names = FALSE, col.names = FALSE)
=======
                sep = "", row.names = FALSE, col.names = FALSE)
>>>>>>> add704465cad8d5bc3d2e2a38b769f5a90d80dac
    x <- eigenstrat(genoFile = "eigenstratG.eg.txt", outFile.Robj = "eigenstrat.result.list",
                    outFile.txt = "eigenstrat.result.txt", rm.marker.index = NULL,
                    rm.subject.index = NULL, miss.val = 9, num.splits = 10,
                    topK = NULL, signt.eigen.level = 0.01, signal.outlier = FALSE,
                    iter.outlier = 5, sigma.thresh = 6)
    x <- x$eigenvectors
    x <- as.matrix(x[,1])
  }
}
snp.tbl <- as.vector(gwas01.Marker[,4])
snp.tbl <- snp.tbl[!is.na(snp.tbl)]
geno.tbl <- as.matrix(geno01[snp.tbl])
geno.tbl[geno.tbl == 0.5] <- 9

cor_mat <- estimate_ss_cor(ref_pcs=tab.pc, ref_genotypes=ref_genotype, link_function='linear')




typeof(geno.tbl)
class(geno.tbl)
geno.tbl


eigenstratG.eg <- matrix(rbinom(3000, 2, 0.5), ncol = 30)
eigenstratG.eg

typeof(eigenstratG.eg)
class(eigenstratG.eg)

write.table(geno.tbl, file = "eigenstratG.eg.txt", quote = FALSE,
            sep = "", row.names = FALSE, col.names = FALSE)

x <- eigenstrat(genoFile = "eigenstratG.eg.txt", outFile.Robj = "eigenstrat.result.list",
                outFile.txt = "eigenstrat.result.txt", rm.marker.index = NULL,
                rm.subject.index = NULL, miss.val = 9, num.splits = 10,
                topK = NULL, signt.eigen.level = 0.01, signal.outlier = FALSE,
                iter.outlier = 5, sigma.thresh = 6)
x <- x$eigenvectors
x <- as.matrix(x[,1], byrow = TRUE)
x

geno.tbl <- as.data.frame(geno.tbl)

y <- estimate_ss_cor(ref_pcs=x, ref_genotypes=geno.tbl, link_function='logit')

typeof(geno.tbl)
class(geno.tbl)

file.remove("eigenstratG.eg.txt", "eigenstrat.result.list", "eigenstrat.result.txt")



system("ls -F")
