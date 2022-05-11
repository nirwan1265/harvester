# Package names
packages <- c("ggplot2", "Rsamtools","GenomicAlignments","rtracklayer","GenomicRanges","AnnotationHub","knitr","gtools","data.table","stringi","GBJ","metap","multtest","Hmisc","devtools","SNPRelate","gdsfmt","dplyr","vcfR","tidyr","AssocTests","SKAT")

# Install packages not yet installed
#installed_packages <- packages %in% rownames(installed.packages())
#if (any(installed_packages == FALSE)) {
#  install.packages(packages[!installed_packages])
#}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

## Read the pathway database:
setwd("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Data for sorghum/sorghum/Pathways/sorghumbicolorcyc")
pathway <- read.delim("pathways.txt", sep ="\t")
pathway <- as.data.frame(t(pathway[,-1]))
colnames(pathway) <- pathway[1,]
pathway <- pathway[-1,]


#Read the combined pvalue files
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Data for sorghum/sorghum/SNP annotation/combined.pvalues")
for(i in sprintf("%02d", c(1,3,7,9,10))){
  assign(paste0("pvalue.chr",i), readRDS(paste0("pvalue.combine.sorghum.chr",i,".RDS")))
}

#Combine all the genes in one file
j <- 1
all.genes = {}
for(i in paste0("pvalue.chr", sprintf("%02d", c(1,3,7,9,10)))){
  d = get(i)
  all.genes <- rbind(all.genes,d)
  assign(i,d)
  j <- j+1
}

all.genes <- as.data.frame(t(all.genes))
all.genes <- all.genes[2,]



#Filtering the pathway genes
x <- as.data.frame(as.matrix(NA))
for(i in 1:ncol(pathway)){
  for(j in 1:nrow(pvalue.chr01)){
    if(rownames(pvalue.chr01)[j] %in% pathway[,i] == TRUE){
      x[j,i] <- rownames(pvalue.chr01)[j] 
    }else {
      x[j,i] <- 0
    }  
  }
}

#Sorting the pathway genes and naming the pathways
sorghum.pathway <- as.data.frame(apply(x,2,sort,decreasing = TRUE))
colnames(sorghum.pathway) <- colnames(pathway)

#Removing empty pathways
sorghum.pathway <- sorghum.pathway[,colSums(sorghum.pathway != 0) > 0]

#Sorting by the highest number of genes
sorghum.pathway <- sorghum.pathway[,order(colSums(sorghum.pathway != 0), decreasing = TRUE)]

write.csv(sorghum.pathway,"sorghum.pathway.csv",  row.names = FALSE, quote = FALSE)




