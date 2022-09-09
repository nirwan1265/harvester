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

