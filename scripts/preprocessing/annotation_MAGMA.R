#Collecting all the genotype SNP names
for(i in sprintf("%02d", 1:10)){
  assign(paste0("x",i), as.data.frame(colnames(get(paste0("geno", i)))))
}


for(i in )
