###PCA analysis required for GBJ
#Need a vcf file format of the hapmap which is converted to GDS format
#TASSEL or PLINK is used for converting hapmap to VCF file format
#Need a directory to  create the gds file. If working on the server, we might need to define this before starting
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Lasky.hapmap/hapmap_vcf/minor_allele_format")

##Reading the vcf files
for(i in sprintf("%02d", 1:10)){
  assign(paste0("vcf.fn",i),paste0("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Lasky.hapmap/raw/vcf/africa_filtered/sb_snpsDryad_sept2013_filter.c",i,".vcf"))
}

##Converting vcf to gds
#A bit time consuming
#Note: for some reason, you cannot run the next step twice if you make an error. you need to delete all this converted gds files, remove all your env variables and do it again.
j <- 1
for(i in paste0("vcf.fn",sprintf("%02d", 1:10))){
  d = get(i)
  snpgdsVCF2GDS(d, paste0("chr",sprintf("%02d",j),".gds"), method = "copy.num.of.ref")
  assign(i,d)
  j = j + 1
}


##Get the GDS file data
for(i in sprintf("%02d", 1:10)){
  assign(paste0("gdsfile",i), snpgdsOpen(paste0("chr",i,".gds")))
}

##LD-based SNP pruning
set.seed(1000)
# Try different LD thresholds for sensitivity analysis but read in a paper somewhere that 0.2 was used for GBJ
j <- 1
for(i in paste0("gdsfile",sprintf("%02d", 1:10))){
  d = get(i)
  assign(paste0("snpset",sprintf("%02d",j)), snpgdsLDpruning(d,ld.threshold = 0.2))
  assign(i,d)
  j = j + 1
}

## Get all selected snp id
j <- 1
for(i in paste0("snpset",sprintf("%02d",1:10))){
  d = get(i)
  assign(paste0("snpset.id",sprintf("%02d",j)), unlist(unname(d)))
  assign(i,d)
  j = j + 1
}

## Run PCA
j <- 1
for(i in paste0("snpset.id",sprintf("%02d",1:10))){
  d = get(i)
  assign(paste0("pca",sprintf("%02d",j)), snpgdsPCA(get(paste0("gdsfile",sprintf("%02d",j))), snp.id = d, num.thread = 10))
  assign(i,d)
  j = j+1
}


#In case there are population information
#https://www.bioconductor.org/packages/devel/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html
#In the case of no prior population information,
#First two shows max variance
#Make a table of eigen values
j <- 1
for(i in paste0("pca",sprintf("%02d",1:10))){
  d = get(i)
  assign(paste0("tab",sprintf("%02d",j)), data.frame(sample.id = d$sample.id,
                                                     EV1 = d$eigenvect[,1],
                                                     EV2 = d$eigenvect[,2],
                                                     EV3 = d$eigenvect[,3],
                                                     EV4 = d$eigenvect[,4],
                                                     EV5 = d$eigenvect[,5],
                                                     EV6 = d$eigenvect[,6],
                                                     EV7 = d$eigenvect[,7],
                                                     EV8 = d$eigenvect[,8],
                                                     EV9 = d$eigenvect[,9],
                                                     EV10 = d$eigenvect[,10],
                                                     stringsAsFactors = FALSE))
  assign(i,d)
  j = j + 1
}

write.csv(tab10, "PCA_LM_sorghum.csv", row.names = F)
# Saving tab
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/R_saved")
# j <- 1
# for(i in paste0("tab",sprintf("%02d", 1:10))){
#   d = get(i)
#   saveRDS(d, paste0("tab",sprintf("%02d",j),".RDS"))
#   j <- j + 1
#   assign(i,d)
# }

system("ls")
tab01 <- readRDS("tab01.RDS")
tab01 <- matrix(as.numeric(tab01), ncol = ncol(tab01))

tab02 <- readRDS("tab02.RDS")
tab02 <- matrix(as.numeric(tab02), ncol = ncol(tab02))

tab03 <- readRDS("tab03.RDS")
tab03 <- matrix(as.numeric(tab03), ncol = ncol(tab03))

tab04 <- readRDS("tab04.RDS")
tab04 <- matrix(as.numeric(tab04), ncol = ncol(tab04))

tab05 <- readRDS("tab05.RDS")
tab05 <- matrix(as.numeric(tab05), ncol = ncol(tab05))

tab06 <- readRDS("tab06.RDS")
tab06 <- matrix(as.numeric(tab06), ncol = ncol(tab06))

tab07 <- readRDS("tab07.RDS")
tab07 <- matrix(as.numeric(tab07), ncol = ncol(tab07))

tab08 <- readRDS("tab08.RDS")
tab08 <- matrix(as.numeric(tab08), ncol = ncol(tab08))

tab09 <- readRDS("tab09.RDS")
tab09 <- matrix(as.numeric(tab09), ncol = ncol(tab09))

tab10 <- readRDS("tab10.RDS")
tab10 <- matrix(as.numeric(tab10), ncol = ncol(tab10))

