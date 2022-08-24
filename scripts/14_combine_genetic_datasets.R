# Loading Data
## ONMI data
omni_uniprot <- read.csv("OMNI_all_filtered.for.uniprot.csv")
head(omni_uniprot)
omni <- omni_uniprot[,c(1,3)]
head(omni)

## RNAseq data
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/RNAseq_lowP")
rnaseq_unirprot <- read.csv("lowP_sorghum_rnaseq_diff.csv")
head(rnaseq_unirprot)
rnaseq <- rnaseq_unirprot[,c(3,8)]
head(rnaseq)
colnames(rnaseq) <- c("pvalue", "Gene")
rnaseq[rnaseq == 0] <- min(rnaseq$pvalue[rnaseq$pvalue > 0])/2

#Combining datasets using Cauchy Combination Test:
omni_rnaseq <- inner_join(omni,rnaseq, by = "Gene")
head(omni_rnaseq)
colnames(omni_rnaseq) <- c("GeneName", "pvalue.omni","pvalue.rnaseq")


## Using function
omni_rnaseq_cct <- cbind(omni_rnaseq, (apply(omni_rnaseq[,c(2,3)],1,acat)))

