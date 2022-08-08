###Loading databases
##Downloading SNP database
#SNP Database 
#http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-52/plants/gff3/sorghum_bicolor/
#Will need to make a folder for all the databases for Sorghum and Maize
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Sorghum.annotation/Phytozome/PhytozomeV12/Sbicolor/annotation")
snp.db <- read.table(file ="Sbicolor_454_v3.1.1.gene_exons.gff3", sep = "\t", header = FALSE)
colnames(snp.db) <- c("Chromosome","Database","Region","Start","End","NA","Strand","NA2","Gene")


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



##SNP and pvalue table
setwd("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/eMAGMA")
chr01.snp <- query.snp.gwas01[,c(1,5)]
write.table(chr01.snp,"chr01.snp.pvalues.txt", sep = " ", row.names = F, col.names = F, quote = F)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


##Loading query SNPs for GWAS from RDS file
##Location in the server: /rsstu/users/r/rrellan/sara/SorghumGEA/results/GLM_20220222 - GLM
##Location in the server: /rsstu/users/r/rrellan/sara/SorghumGEA/results/GLM_20220224 - GLM
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/GWAS.results/MLM")
##Reading the gwas RDS files
for(i in sprintf("%02d", 1:10)){
  assign(paste0("query.snp.gwas", i) , readRDS(file = paste0("mlm_sol_VL_",i,"_20220307_19_13.RDS")))
}

##Subsetting the required columns for analysis
for(i in sprintf("%02d", 1:10)){
  assign(paste0("query.snp.gwas",i), data.frame(get(paste0("query.snp.gwas",i))[3]))
  assign(paste0("query.snp.gwas",i), get(paste0("query.snp.gwas",i))[,c(2,3,4,6,7)])
}

##Sub-setting chromosomes from gene annotation database
for(i in sprintf("%02d", 1:10)){
  assign(paste0("x",i), paste0("Chr",i))
}
for(i in sprintf("%02d", 1:10)){
  assign(paste0("db.",i), snp.db[which(snp.db$Chromosome == get(paste0("x",i))), ])
}


##Adding z value column
for(i in paste0("gwas", sprintf("%02d", 1:10))){
  d = get(i)
  d$zstat = as.data.frame(apply(d,1,zval))
  assign(i,d)
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


###Annotation the SNPs
##Making GRanges for Database 
for(i in sprintf("%02d", 1:10)){
  assign(paste0("gr.db", i) , GRanges(seqnames = paste0("chr",i), ranges = IRanges(start = get(paste0("db.",i))[,"Start"], end = get(paste0("db.",i))[,"End"]), strand = get(paste0("db.",i))[,"Strand"], Region = get(paste0("db.",i))[,"Region"], Gene = get(paste0("db.",i))[,"Gene"]))
}



##Making GRanges for gwas Query
#Changing the position column to numeric and removing the first row
for(i in paste0("query.snp.gwas", sprintf("%02d", 1:10))){
  d = get(i)
  d = d[-1,]
  d$Pos = as.numeric(d$MLM_Stats.Pos)
  assign(i,d)
}

for(i in sprintf("%02d", 1:10)){
  assign(paste0("gr.q", i) , GRanges(seqnames = paste0("chr",i), ranges = IRanges(start = get(paste0("query.snp.gwas",i))[,"Pos"], width = 1, fstat = get(paste0("query.snp.gwas",i))[,"MLM_Stats.F"], Marker = get(paste0("query.snp.gwas",i))[,"MLM_Stats.Marker"],pvalue = get(paste0("query.snp.gwas",i))[,"MLM_Stats.p"])))
}


#Finding the Overlaps
for(i in sprintf("%02d", 1:10)){
  assign(paste0("common",i), as.data.frame(findOverlapPairs(get(paste0("gr.db",i)), get(paste0("gr.q",i)))))
}


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


###Data prep
##Filtering out the gene, fstat, Marker and pvalue columns
##Making a new dummy table
for(i in sprintf("%02d", 1:10)){
  assign(paste0("gwas",i), get(paste0("common",i)))
}

# #Saving RDS
# j <- 1
# for(i in paste0("common", sprintf("%02d", 1:10))){
#   d = get(i)
#   saveRDS(d, paste0("common",sprintf("%02d" , j),".RDS"))
#   assign(i,d)
#   j <- j+1
# }



##Filter table having only gene
for(i in paste0("gwas", sprintf("%02d", 1:10))){
  d=get(i)
  d <- d[which(d$first.X.Region == "gene"), ]
  assign(i,d)
}

##Filtering out the required columns for analysis
for(i in sprintf("%02d", 1:10)){
  assign(paste0("gwas",i), data.frame(get(paste0("gwas",i)))[,c(7,15,16,17)])
}

##Renaming columns
for(i in paste0("gwas", sprintf("%02d", 1:10))){
  d=get(i)
  colnames(d) = c("Gene","fstat","Marker","pvalue")
  assign(i,d)
}

##Sorting a/c gene name
#Not required, but might be useful for other database or when combining all results
# for(i in paste0("gwas", sprintf("%02d", 1:10))){
#   d=get(i)
#   d <- d[gtools::mixedorder(d$Gene), ]
#   assign(i,d)
# }

##Table with Z-stat values
for(i in paste0("gwas",sprintf("%02d",1:10))){
  d=get(i)
  assign(paste0(i,".zstat"), dcast(setDT(d), Gene~rowid(Gene, prefix = "zstat"), value.var = "zstat"))
  assign(i,d)
}
#----> continue
##Adding gene names
#Table with all the gene names from all the chromosome
#Can use this table for all the other tables with Marker and pvalue data
for(i in paste0("gwas", sprintf("%02d", 1:10))){
  d=get(i)
  assign(paste0(i,".gene.names"), get(paste0(i,".zstat"))[,1])
  assign(i,d)
}

for(i in paste0("gwas", sprintf("%02d", 1:10),".gene.names")){
  d=get(i)
  d <-  as.data.frame(apply(d,1,split.names))
  assign(i,d)
}

#----> continue
#Adding gene names
j <- 1
for(i in paste0("gwas",sprintf("%02d",1:10),".zstat")){
  d <- get(i)
  d[,1] <- get(paste0("gwas",sprintf("%02d",j),".gene.names"))[,1]
  assign(i,d)
  j = j + 1
}
#Data ordering
for(i in paste0("gwas",sprintf("%02d", 1:10),".zstat")){
  d <- get(i)
  d <- as.data.frame(t(d))
  x <- d[1,]
  d <- d[-1,]
  d <- apply(d,2,sort)
  d <- as.data.frame(stri_list2matrix(d, byrow = FALSE))
  colnames(d) <- x
  assign(i,d)
}


##Table with SNP markers
for(i in paste0("gwas",sprintf("%02d",1:10))){
  d=get(i)
  assign(paste0(i,".Marker"), dcast(setDT(d), Gene~rowid(Gene, prefix = "Marker"), value.var = "Marker"))
  assign(i,d)
}
#Adding gene names
j <- 1
for(i in paste0("gwas",sprintf("%02d",1:10),".Marker")){
  d <- get(i)
  d[,1] <- get(paste0("gwas",sprintf("%02d",j),".gene.names"))[,1]
  assign(i,d)
  j = j + 1
}
#Data ordering
for(i in paste0("gwas",sprintf("%02d", 1:10),".Marker")){
  d <- get(i)
  d <- as.data.frame(t(d))
  x <- d[1,]
  d <- d[-1,]
  d <- apply(d,2,sort)
  d <- as.data.frame(stri_list2matrix(d, byrow = FALSE))
  colnames(d) <- x
  assign(i,d)
}

##Table with pvalue values
for(i in paste0("gwas",sprintf("%02d",1:10))){
  d=get(i)
  assign(paste0(i,".pvalue"), dcast(setDT(d), Gene~rowid(Gene, prefix = "pvalue"), value.var = "pvalue"))
  assign(i,d)
}
#Adding gene names
j <- 1
for(i in paste0("gwas",sprintf("%02d",1:10),".pvalue")){
  d <- get(i)
  d[,1] <- get(paste0("gwas",sprintf("%02d",j),".gene.names"))[,1]
  assign(i,d)
  j = j + 1
}
#Data ordering
for(i in paste0("gwas",sprintf("%02d", 1:10),".pvalue")){
  d <- get(i)
  d <- as.data.frame(t(d))
  x <- d[1,]
  d <- d[-1,]
  d <- apply(d,2,sort)
  d <- as.data.frame(stri_list2matrix(d, byrow = FALSE))
  colnames(d) <- x
  assign(i,d)
}


#saving objects
# j <- 1
# for(i in paste0("gwas", sprintf("%02d", 1:10),".fstat")){
#   d = get(i)
#   saveRDS(d, paste0("gwas",sprintf("%02d" , j),".fstat.RDS"))
#   assign(i,d)
#   j <- j+1
# }
# j <- 1
# for(i in paste0("gwas", sprintf("%02d", 1:10),".gene.names")){
#   d = get(i)
#   saveRDS(d, paste0("gwas",sprintf("%02d" , j),".gene.names.RDS"))
#   assign(i,d)
#   j <- j+1
# }
# j <- 1
# for(i in paste0("gwas", sprintf("%02d", 1:10),".Marker")){
#   d = get(i)
#   saveRDS(d, paste0("gwas",sprintf("%02d" , j),".Marker.RDS"))
#   assign(i,d)
#   j <- j+1
# }
# j <- 1
# for(i in paste0("gwas", sprintf("%02d", 1:10),".pvalue")){
#   d = get(i)
#   saveRDS(d, paste0("gwas",sprintf("%02d" , j),".pvalue.RDS"))
#   assign(i,d)
#   j <- j+1
# }
# 


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


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


###Pvalue combination
##Pre-processing for p-value combination
for(i in sprintf("%02d", 1:10)){
  assign(paste0("tab.pc",i), get(paste0("tab",i))[,c(2:3)])
}

## Loading Numerical hapmap genotype file
#Numerical file should contain the only the taxa present in the phenotype.
#this is done by filtering with taxa in tassel
#Put the taxa name on the bar in tassel and filter

#Save as numerical in TASSEL
#IMPORTANT: This can be done in TASSEL. Remove <marker> and <numerical> text from the transposed txt file before loading 
#pre processing step
#Remove the first Element <Marker>
# $ cut -f2- chr1.txt > chr1.num.txt
#Remove the first row 
# $ sed -e '1d' < chr1.num.txt > chr1.1num.txt
#Separate the file in header and body
#awk 'FNR <=1' sb_snpsDryad_sept2013_filter.c01_numerical.txt > header1.txt
#sed '1d' sb_snpsDryad_sept2013_filter.c01_numerical.txt > body1.txt
#Script to change NA to 9, 0 to 2, 1 to 0 and 0.5 to 1 and then 2.5 to 1 because it replaces 0.5 as 2.5 
# $ perl -pi -e 's/NA/9/g' body1.txt
#combine the header and body
# cat header1.txt body1.txt > originalfile.txt


##Reading the genotype files
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Lasky.hapmap/hapmap_numerical/tassel_format")

for(i in sprintf("%02d", 1)){
  assign(paste0("geno",i), read.table(file = paste0("sb_snpsDryad_sept2013_filter.c",i,"_numerical.txt"), header = TRUE, sep = "\t"))
}


##Saving the genotype file as  RDS
# j <- 1
# for(i in paste0("geno", sprintf("%02d", 1:10))){
#   d = get(i)
#   saveRDS(d, paste0("geno",sprintf("%02d" , j),".RDS"))
#   assign(i,d)
#   j <- j+1
# }

#save.image(file="sorghum_VL_omnibus.RData")
#load(file="sorghum_VL_omnibus.RData")


##Combination tests using GBJ package

#Convert tab.pcs to matrix
for(i in paste0("tab", sprintf("%02d", 1:10))){
  d = get(i)
  d <- as.matrix(d[,-1])
  assign(i,d)
}


#Convert zstat dataframe to numeric
for(i in paste0("gwas", sprintf("%02d", 1:10), ".zstat")){
  d <- get(i)
  d <- as.data.frame(lapply(d, as.numeric))
  assign(i,d)
}

##Saving the tab.pcs file as  RDS
# j <- 1
# for(i in paste0("gwas", sprintf("%02d", 1:10))){
#   d = get(i)
#   saveRDS(d, paste0("tab.pc",sprintf("%02d" , j),".RDS"))
#   assign(i,d)
#   j <- j+1
# }



#Using Global Function
# for(i in sprintf("%02d", 1:10)){
#   assign(paste0("pvalue.combine",i), pvalue.combine(get(paste0("gwas",i,".fstat")), get(paste0("gwas",i,".Marker")), get(paste0("gwas",i,".pvalue")), get(paste0("geno",i)), get(paste0("tab.pc",i))))
# }
# 
# 
# ##Adding gene(row) and test(column) names
# j <- 10
# for(i in paste0("pvalue.combine",sprintf("%02d", 10))){
#   d = get(i)
#   row.names(d) <- get(paste0("gwas", sprintf("%02d", j), ".gene.names"))[,1]
#   colnames(d) <- c("GBJ","GHC","minP","OMNI")
#   assign(i,d)
#   j = j + 1
# }









#number of phenotype should match genotype file. Columns should match. 

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#Phenotype file for SKAT
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Phenotype")
pheno <- read.table("VL_P.txt")
pheno_name <- rownames(pheno) 
pheno_name <- pheno_name[-c(1,2)]
pheno <- as.matrix(as.integer(pheno[-c(1,2),]))
rownames(pheno) <- pheno_name

#Genotype file for SKAT
#filter by taxa from TASSEL, i.e., only the taxa which are used in the phenotype should be present in the genotype file 
#Do the same as before


##Reading the genotype files
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Lasky.hapmap/genotype_filtered_by_phenotype/minor_allele_format")
#edit the files like the one from geno
for(i in sprintf("%02d", 1:10)){
  assign(paste0("geno",i), read.table(file = paste0("geno_f_pheno",i,".txt"), header = TRUE, sep = "\t"))
}


##Saving the genotype file as  RDS
# j <- 1
# for(i in paste0("geno", sprintf("%02d", 1:10))){
#   d = get(i)
#   saveRDS(d, paste0("geno_f",sprintf("%02d" , j),".RDS"))
#   assign(i,d)
#   j <- j+1
# }





#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


#Running the Analysis:
#Using the pvalue.combination function:
for(i in sprintf("%02d", 6:10)){
  assign(paste0("pvalue.combine",i), pvalue.combine(get(paste0("gwas",i,".zstat")), get(paste0("gwas",i,".Marker")), get(paste0("gwas",i,".pvalue")), get(paste0("geno",i)), get(paste0("tab",i))))
}

#Adding names
j = 1
for(i in paste0("pvalue.combine",sprintf("%02d", 1:10))){
  d = get(i)
  a <- get(paste0("gwas",sprintf("%02d", j), ".gene.names"))
  a <- unlist(a)
  row.names(d) <- a
  j = j + 1
  assign(i,d)
}


#Column names
for(i in paste0("pvalue.combine",sprintf("%02d", 1:10))){
  d = get(i)
  colnames(d) <- c("GBJ","GHC","minP","SKAT","OMNI_ItoIV","CCT_ItoIV")
  assign(i,d)
}

#SAVING the files
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/pvalues.combination")
# j = 1
# for(i in paste0("pvalue.combine",sprintf("%02d", 1:10))){
#   d = get(i)
#   write.csv(d, paste0("pvalue.combine",sprintf("%02d", j),".csv"))
#   j = j + 1
#   assign(i,d)
# }


#Reading the files
# for(i in sprintf("%02d", 1:5)){
#   assign(paste0("pvalue.combine",i), read.csv(paste0("pvalue.combine",i,".csv")))
# }

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


##Plotting the number of SNPs for each gene
par(mfrow = c(5, 5))

for(i in paste0("gwas",sprintf("%02d", 1:10),".Marker")){
  d = get(i)
  elements <- as.data.frame(colSums(!is.na(d)))
  colnames(elements) <- "Number.of.SNPs"
  hist(elements$Number.of.SNPs, main = "Distribution of SNPs",
       xlab = "Number of SNPs", 
       col = rainbow(14),
       breaks = max(elements), #highest SNPs
       ylim = c(1,1000),
       labels = TRUE
  )
  assign(i,d)
}


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


#Filtering only the significant genes
filtered_genes_OMNI <- NULL
filtered_genes_CCT <- NULL
filtered_genes_GBJ <- NULL
for(i in paste0("pvalue.combine",sprintf("%02d", 1:10))){
  d = get(i)
  x <- dplyr::filter(d, OMNI_ItoIV < 0.05)
  y <- dplyr::filter(d, CCT_ItoIV < 0.05)
  z <- dplyr::filter(d, GBJ < 0.05)
  filtered_genes_OMNI <- rbind(filtered_genes_OMNI,x)
  filtered_genes_CCT <- rbind(filtered_genes_CCT,y)
  filtered_genes_GBJ <- rbind(filtered_genes_GBJ,z)
  assign(i,d)
}

