###Loading databases
##Downloading SNP database
#SNP Database 
#http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-52/plants/gff3/sorghum_bicolor/
#Will need to make a folder for all the databases for Sorghum and Maize
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Sorghum.annotation/Phytozome/PhytozomeV12/Sbicolor/annotation")
snp.db <- read.table(file ="Sbicolor_454_v3.1.1.gene_exons.gff3", sep = "\t", header = FALSE)
colnames(snp.db) <- c("Chromosome","Database","Region","Start","End","NA","Strand","NA2","Gene")



##SNP and pvalue table
setwd("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/eMAGMA")
chr01.snp <- query.snp.gwas01[,c(1,5)]
write.table(chr01.snp,"chr01.snp.pvalues.txt", sep = " ", row.names = F, col.names = F, quote = F)


##Sub-setting chromosomes from gene annotation database
for(i in sprintf("%02d", 1:10)){
  assign(paste0("x",i), paste0("Chr",i))
}
for(i in sprintf("%02d", 1:10)){
  assign(paste0("db.",i), snp.db[which(snp.db$Chromosome == get(paste0("x",i))), ])
}

