## Read the pathway database:
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Data for sorghum/sorghum/Pathways/sorghumbicolorcyc/7.0/data")
pathway <- read.delim("pathways.txt", sep ="\t")
pathway <- as.data.frame(t(pathway[,-1]))
colnames(pathway) <- pathway[1,]
pathway <- pathway[-1,]

#Read the combined pvalue files
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Data for sorghum/sorghum/SNP annotation/combined.pvalues")
for(i in sprintf("%02d", 3)){
  assign(paste0("pvalue.chr",i), readRDS(paste0("pvalue.combine.sorghum.chr",i,".RDS")))
}

pvalue.chr01 <- as.data.frame(t(pvalue.chr03))
pvalue.chr01 <- pvalue.chr01[2,]
pvalue.chr01


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
rownames(sorghum.pathway) <- colnames(pathway)

#Removing empty pathways
sorghum.pathway <- sorghum.pathway[,colSums(sorghum.pathway != 0) > 0]

#Sorting by the highest number of genes
sorghum.pathway <- sorghum.pathway[,order(colSums(sorghum.pathway != 0), decreasing = TRUE)]

write.csv(sorghum.pathway,"sorghum.pathway.csv",  row.names = FALSE, quote = FALSE)




