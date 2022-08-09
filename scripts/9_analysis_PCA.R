###PCA analysis required for GBJ
#Need a vcf file format of the hapmap which is converted to GDS format
#TASSEL or PLINK is used for converting hapmap to VCF file format
#Need a directory to  create the gds file. If working on the server, we might need to define this before starting
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Lasky.hapmap/hapmap_vcf/minor_allele_format")

##Reading the vcf files
for(i in sprintf("%02d", 1:10)){
  assign(paste0("vcf.fn",i),paste0("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Lasky.hapmap/hapmap_vcf/minor_allele_format/sb_snpsDryad_sept2013_filter.c",i,".vcf"))
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
  assign(paste0("pca",sprintf("%02d",j)), snpgdsPCA(get(paste0("gdsfile",sprintf("%02d",j))), snp.id = d, num.thread = 2))
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
                                                     stringsAsFactors = FALSE))
  assign(i,d)
  j = j + 1
}


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

