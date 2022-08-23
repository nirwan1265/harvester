#Running the Analysis:
#Using the pvalue.combination function:
for(i in sprintf("%02d", 1:10)){
  assign(paste0("pvalue.combine",i), pvalue.combine(get(paste0("gwas",i,".zstat")), get(paste0("gwas",i,".Marker")), get(paste0("gwas",i,".pvalue")), get(paste0("geno",i)), get(paste0("tab",i))))
}

# Saving as R objects
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/R_saved")
# saveRDS(pvalue.combine01,"pvalue.combine01.RDS")
# saveRDS(pvalue.combine02,"pvalue.combine02.RDS")
# saveRDS(pvalue.combine03,"pvalue.combine03.RDS")
# saveRDS(pvalue.combine04,"pvalue.combine04.RDS")
# saveRDS(pvalue.combine05,"pvalue.combine05.RDS")
# saveRDS(pvalue.combine06,"pvalue.combine06.RDS")
# saveRDS(pvalue.combine07,"pvalue.combine07.RDS")
# saveRDS(pvalue.combine08,"pvalue.combine08.RDS")
# saveRDS(pvalue.combine09,"pvalue.combine09.RDS")
# saveRDS(pvalue.combine10,"pvalue.combine10.RDS")


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

