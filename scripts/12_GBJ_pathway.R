
# Read the pathway database:
# To convert to ENSEMBL or NCBI naming system, change Sobic. to SORBI_3
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Sorghum_pathway_database/pythozome/sorghumbicolorcyc/7.0/data")
pathway <- read.csv("pathways_ensembl.csv")
pathway <- as.data.frame(t(pathway[,-1]))
colnames(pathway) <- pathway[1,]
pathway <- pathway[-1,]



#Read the combined pvalue files
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/pvalues.combination")
# for(i in sprintf("%02d", 1:10)){
#   assign(paste0("pvalue.combine",i), readRDS(paste0("pvalue.combine",i,".RDS")))
# }
 
#All genes:
# write.csv(GBJ_all,"GBJ_all.csv")
# write.csv(GHC_all,"GHC_all.csv")
# write.csv(minP_all,"minP_all.csv")
# write.csv(SKAT_all,"SKAT_all.csv")
# write.csv(OMNI_all,"OMNI_all.csv")
# write.csv(CCT_all,"CCT_all.csv")
# Read files


#Filtering the pathway genes
x <- as.data.frame(as.matrix(NA))
for(i in 1:ncol(pathway)){
  for(j in 1:nrow(OMNI_all)){
    if(rownames(filtered_genes_OMNI)[j] %in% pathway[,i] == TRUE){
      x[j,i] <- rownames(filtered_genes_OMNI)[j] 
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

#SAVING AND READING pathway files
write.csv(sorghum.pathway,"pathway_OMNI.csv")


####################################################################################
####################################################################################

# pvalue combination for the pathways






