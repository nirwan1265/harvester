library(ggvenn)
library(VennDiagram)
library(venn)

# Recover data
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/pvalues.combination")



##Plotting the number of SNPs for each gene
# par(mfrow = c(2, 5))
# 
# for(i in paste0("gwas",sprintf("%02d", 1:10),".Marker")){
#   d = get(i)
#   elements <- as.data.frame(colSums(!is.na(d)))
#   colnames(elements) <- "Number.of.SNPs"
#   hist(elements$Number.of.SNPs, main = "Distribution of SNPs",
#        xlab = "Number of SNPs", 
#        col = rainbow(14),
#        breaks = max(elements), #highest SNPs
#        ylim = c(1,1000),
#        labels = TRUE
#   )
#   assign(i,d)
# }
# 

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Read the files
for(i in sprintf("%02d", 1:10)){
  assign(paste0("pvalue.combine",i), read.csv(paste0("pvalue.combine",i,".csv")))
}
j <- 1
for(i in paste0("pvalue.combine",sprintf("%02d", 1:10))){
  d = get(i)
  colnames(d) <- c("Gene","GBJ","GHC","minP","SKAT","OMNI_ItoIV","CCT_ItoIV")
  assign(i,d)
}

#Filtering ALL GENES
OMNI_all <- NULL
CCT_all <- NULL
GBJ_all <- NULL
SKAT_all <- NULL
GHC_all <- NULL
minP_all <- NULL
for(i in paste0("pvalue.combine",sprintf("%02d", 1:10))){
  d = get(i)
  u <- dplyr::filter(d, GBJ < 5)
  v <- dplyr::filter(d, GHC < 5)
  w <- dplyr::filter(d, minP < 5)
  x <- dplyr::filter(d, SKAT < 5)
  y <- dplyr::filter(d, OMNI_ItoIV < 5)
  z <- dplyr::filter(d, CCT_ItoIV < 5)
  GBJ_all <- rbind(GBJ_all,u)
  GHC_all <- rbind(GHC_all,v)
  minP_all <- rbind(minP_all,w)
  SKAT_all <- rbind(SKAT_all,x)
  OMNI_all <- rbind(OMNI_all,y)
  CCT_all <- rbind(CCT_all,z)
  assign(i,d)
}
x <- c("GBJ","GHC","minP","SKAT","OMNI_ItoIV","CCT_ItoIV")
y <- c("GBJ_all","GHC_all","minP_all","SKAT_all","OMNI_all","CCT_all")



#Filtering only the significant genes
filtered_genes_OMNI <- NULL
filtered_genes_CCT <- NULL
filtered_genes_GBJ <- NULL
filtered_genes_SKAT <- NULL
filtered_genes_GHC <- NULL
filtered_genes_minP <- NULL
for(i in paste0("pvalue.combine",sprintf("%02d", 1:10))){
  d = get(i)
  u <- dplyr::filter(d, GBJ < 0.05)
  v <- dplyr::filter(d, GHC < 0.05)
  w <- dplyr::filter(d, minP < 0.05)
  x <- dplyr::filter(d, SKAT < 0.05)
  y <- dplyr::filter(d, OMNI_ItoIV < 0.05)
  z <- dplyr::filter(d, CCT_ItoIV < 0.05)
  filtered_genes_GBJ <- rbind(filtered_genes_GBJ,u)
  filtered_genes_GHC <- rbind(filtered_genes_GHC,v)
  filtered_genes_minP <- rbind(filtered_genes_minP,w)
  filtered_genes_SKAT <- rbind(filtered_genes_SKAT,x)
  filtered_genes_OMNI <- rbind(filtered_genes_OMNI,y)
  filtered_genes_CCT <- rbind(filtered_genes_CCT,z)
  assign(i,d)
}

## Save the files
write.csv(filtered_genes_CCT,"filtered_genes_CCT.csv", row.names = FALSE)
write.csv(filtered_genes_GBJ,"filtered_genes_GBJ.csv", row.names = FALSE)
write.csv(filtered_genes_GHC,"filtered_genes_GHC.csv", row.names = FALSE)
write.csv(filtered_genes_minP,"filtered_genes_minP.csv", row.names = FALSE)
write.csv(filtered_genes_OMNI,"filtered_genes_OMNI.csv", row.names = FALSE)
write.csv(filtered_genes_SKAT,"filtered_genes_SKAT.csv", row.names = FALSE)

system("pwd")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#Filtered genes:
filtered_genes_OMNI <- read.csv("filtered_genes_OMNI.csv")
filtered_genes_CCT <- read.csv("filtered_genes_CCT.csv")
filtered_genes_GBJ <- read.csv("filtered_genes_GBJ.csv")
filtered_genes_SKAT <- read.csv("filtered_genes_SKAT.csv")
filtered_genes_minP <- read.csv("filtered_genes_minP.csv")
filtered_genes_GHC <- read.csv("filtered_genes_GHC.csv")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


# Venn diagram
x <- list(
  GBJ = as.vector(filtered_genes_GBJ[,1]),
  GHC = as.vector(filtered_genes_GHC[,1]),
  minP = as.vector(filtered_genes_minP[,1]),
  SKAT = as.vector(filtered_genes_SKAT[,1]),
  OMNI <- as.vector(filtered_genes_OMNI[,1]),
  CCT <- as.vector(filtered_genes_CCT[,1])
)
y <- list(
  OMNI <- as.vector(filtered_genes_OMNI[,1]),
  CCT <- as.vector(filtered_genes_CCT[,1])
)
z <- list(
  GBJ = as.vector(filtered_genes_GBJ[,1]),
  GHC = as.vector(filtered_genes_GHC[,1]),
  minP = as.vector(filtered_genes_minP[,1]),
  SKAT = as.vector(filtered_genes_SKAT[,1])
)
venn(x, snames = c("GBJ","GHC","minP","SKAT","OMNI","CCT"), ilabels = TRUE, ellipse = TRUE, zcolor = "style", opacity = 0.2, plotsize = 15, borders = TRUE, box = T, ggplot = TRUE)
venn(y, snames = c("OMNI","CCT"), ilabels = TRUE, ellipse = TRUE, zcolor = "style", opacity = 0.2, plotsize = 15, borders = TRUE, box = T, ggplot = TRUE)
venn(z, snames = c("GBJ","GHC","minP","SKAT"), ilabels = TRUE, ellipse = TRUE, zcolor = "style", opacity = 0.2, plotsize = 15, borders = TRUE, box = T, ggplot = TRUE)



