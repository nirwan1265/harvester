###Loading databases
##Downloading SNP database
#SNP Database 
#http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-52/plants/gff3/sorghum_bicolor/
#Will need to make a folder for all the databases for Sorghum and Maize
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Sorghum.annotation/ensemblgenomes")


library(ape)
x <- read.gff("Sorghum_bicolor.Sorghum_bicolor_NCBIv3.54.chromosome.1.gff3",na.strings = c(".", "?"), GFF3 = TRUE)
x <- tr2g_gff3("Sorghum_bicolor.Sorghum_bicolor_NCBIv3.54.chromosome.1.gff3", write_tr2g = FALSE, get_transcriptome = FALSE, save_filtered_gff = FALSE)

BiocManager::install("BUSpaRse")
library(BUSpaRse)
##Sub-setting chromosomes from gene annotation database

# For NCBI genes
# setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Sorghum.annotation/Phytozome/PhytozomeV12/Sbicolor/annotation")
# snp.db_2 <- read.table(file ="Sbicolor_454_v3.1.1.gene.gff3", sep = "\t", header = TRUE)
# colnames(snp.db_2) <- c("Chromosome","Database","Region","Start","End","NA","Strand","NA2","Gene")
# for(i in 1:10){
#   assign(paste0("x",i), paste0(i))
# }
# for(i in sprintf("%02d", 1:10)){
#   assign(paste0("db.",i,"_2"), snp.db_2[which(snp.db_2$Chromosome == get(paste0("x",i))), ])
# }
# 

j <- 1
for(i in 1:10 ){
  assign(paste0("db.",sprintf("%02d", j)), read.table(file = paste0("Sorghum_bicolor.Sorghum_bicolor_NCBIv3.54.chromosome.",i,".gff3"), sep = "\t", header = TRUE))
  j <-  j + 1
}

for(i in paste0("db.", sprintf("%02d", 1:10))){
  d = get(i)
  colnames(d) <- c("Chromosome","Database","Region","Start","End","NA","Strand","NA2","Gene")
  assign(i,d)
}



# For maize
setwd("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Maize/Maize.annotation")
system("ls")
maize_ref <- read.table("Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3", sep ="\t", header = T)

chr <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10")
for (i in chr){
  assign(paste0("ref_",i), maize_ref[which(maize_ref[1] == i), ])
}

for(i in paste0("ref_chr", 1:10)){
  d = get(i)
  d = d[which(d[3] == "gene"), ]
  colnames(d) <- c("Chromosome","Database","Region","Start","End","NA","Strand","NA2","Gene")
  assign(i,d)
}
