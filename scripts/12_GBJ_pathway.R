
# Read the pathway database:
# To convert to ENSEMBL or NCBI naming system, change Sobic. to SORBI_3
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Sorghum_pathway_database/pythozome/sorghumbicolorcyc/7.0/data")
pathway <- read.table("pathways_ensembl.txt", sep = "\t", fill = TRUE)

pathway <- as.data.frame(t(pathway[,-1]))
colnames(pathway) <- pathway[1,]
pathway <- pathway[-1,]
rownames(pathway) <- pathway[,1]
pathway <- pathway[,-1]


#Read the combined pvalue files
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/pvalues.combination")
# for(i in sprintf("%02d", c(1,3,7,9,10))){
#   assign(paste0("pvalue.chr",i), readRDS(paste0("pvalue.combine.sorghum.chr",i,".RDS")))
# }

filtered_genes_OMNI <- pvalue.combine01_omni
rownames(filtered_genes_OMNI) <- filtered_genes_OMNI[,2]

ncol(pathway)
nrow(filtered_genes_OMNI)

#Filtering the pathway genes
x <- as.data.frame(as.matrix(NA))
for(i in 1:ncol(pathway)){
  for(j in 1:nrow(filtered_genes_OMNI)){
    if(rownames(filtered_genes_OMNI)[j] %in% pathway[,i] == TRUE){
      x[j,i] <- rownames(filtered_genes_OMNI)[j] 
    }else {
      x[j,i] <- 0
    }  
  }
}

rownames(filtered_genes_OMNI)
[j] %in% pathway[,i]

#Sorting the pathway genes and naming the pathways
sorghum.pathway <- as.data.frame(apply(x,2,sort,decreasing = TRUE))
colnames(sorghum.pathway) <- colnames(pathway)

#Removing empty pathways
sorghum.pathway <- sorghum.pathway[,colSums(sorghum.pathway != 0) > 0]

#Sorting by the highest number of genes
sorghum.pathway <- sorghum.pathway[,order(colSums(sorghum.pathway != 0), decreasing = TRUE)]
head(sorghum.pathway)

write.csv(sorghum.pathway,"pathway.csv")



