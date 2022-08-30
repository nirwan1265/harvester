# Load MAGMA data:
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/pvalues.combination/MAGMA")
magma01 <- read.table("chr01.magma.multi.snpwise.meantop.linreg.txt.genes.out", sep = "", header = TRUE)
magma01_p.multi <- as.data.frame(magma01$P_MULTI)
magma01_p.multi$GENE <- magma01$GENE
row.names(magma01_p.multi) <- magma01$GENE
magma01_snpmean <- as.data.frame(magma01$P_SNPWISE_MEAN)

# Load p-values:
pvalue.combine01

pvalue.combine01_omni <- as.data.frame(pvalue.combine01$OMNI_ItoIV)
pvalue.combine01_omni$GENE <- row.names(pvalue.combine01)


#Combine the dataframes
omni.magma.multi <- inner_join(pvalue.combine01_omni,magma01_p.multi)

