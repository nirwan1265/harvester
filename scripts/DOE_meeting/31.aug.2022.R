# Reading GBJ package files
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/pvalues.combination")
omni_all <- read.csv("OMNI_all.csv")
omni_all <- omni_all[,c(1,6)]
gbj_all <- read.csv("GBJ_all.csv")
gbj_all <- gbj_all[,c(1,2)]
ghc_all <- read.csv("GHC_all.csv")
ghc_all <- ghc_all[,c(1,3)]
minP_all <- read.csv("minP_all.csv")
minP_all <- minP_all[,c(1,4)]
skat_all <- read.csv("SKAT_all.csv")
skat_all <- skat_all[,c(1,5)]
cct_all <- read.csv("cct_all.csv")
cct_all <- cct_all[,c(1,7)]


#Reading Pathways
pathway <- read.csv("pathway_OMNI.csv")

#Reading eMAMGA files
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Github/eMAGMA")
magma_mean_top.1 <- read.delim("pval.comb.txt.genes.out", sep="")
magma_mean_top.1 <- magma_mean_top.1[,c(1,9,10,11)]
colnames(magma_mean_top.1)[1] <- "X"
magma_top <- read.delim("pval.comb2.txt.genes.out", sep="")
magma_top <- magma_top[,c(1,9)]
colnames(magma_top)[1] <- "X"




# All common genes across the analysis
all_list <- list(gbj_all, ghc_all,minP_all, skat_all, omni_all,cct_all, magma_mean_top.1, magma_top)
multi_inner <- Reduce(
  function(x, y ) merge(x, y), 
  all_list
)
colnames(multi_inner) <- c("GENE","GBJ[1]","GHC[2]","minP[3]","SKAT[4]","OMNI(1_to_4)","Cauchy(1_to_4)","magma_combine(5_to_6)","SNPWISE_mean[5]","SNPWISE_Top,1[6]","SNPWISE_Top[7]")
multi_inner <- multi_inner[,c(1,2,3,4,5,9,10,11,6,7,8)]
row.names(multi_inner) <- multi_inner[,1]
multi_inner <- multi_inner[,-1]

#Doing acat for the first 5 dataset:
cauchy_1_to_5 <- as.data.frame(apply(multi_inner,1,acat, n = 5))

#Adding to the multi_inner
multi_inner$Cauchy_1_to_5 <- cauchy_1_to_5[,1]
colnames(multi_inner)[11] <- "Cauchy(1_to5)"
multi_inner <- multi_inner[,c(1,2,3,4,5,6,7,8,9,11,10)]


#Removing row if they have only one SNP

keep <- apply(multi_inner[1:4], 1, function(x) length(unique(x[!is.na(x)])) != 1)
multi_inner_unique <- multi_inner[keep, ]



#Saving the whole file
write.csv(multi_inner,"all_pval.combination.csv")
write.csv(multi_inner_unique,"all.unique_pval.combination.csv")
