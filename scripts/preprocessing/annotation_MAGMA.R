################################################################################
#SNP FILE
################################################################################
#query.snp has it all
j <- 1
for(i in paste0("query.snp.gwas", sprintf("%02d", 1:10))){
  d = get(i)
  assign(paste0("x", sprintf("%02d", j)), d[,c(1,2,3)])
  j = j + 1
  assign(i,d)
}

#Combining all the base pairs
snp.file <- NULL
for(i in paste0("x", sprintf("%02d", 1:10))){
  d = get(i)
  snp.file = rbind(snp.file,d)
  assign(i,d)
}


################################################################################
#GENE FILE
################################################################################
# gene ID, chromosome, start site, stop site
x <- NULL
j <-  1
for(i in paste0("db.", sprintf("%02d", 1:10))){
  d = get(i)
  assign(paste0("x",sprintf("%02d", j)), d[, c(9,4,5,1,3)])
  j = j + 1
  assign(i,d)
}

for(i in paste0("x", sprintf("%02d", 1:10))){
  d = get(i)
  d = d[which(d$Region == "gene"), ]
  assign(i,d)
}


#Spliting gene names
split.names <- function(x,split){
  split.genename <- unlist(strsplit(x, split = ';', fixed = TRUE))[2]
  split.genename2 <- unlist(strsplit(split.genename, split = "=", fixed = TRUE))[2]
  return(split.genename2)
}


for(i in paste0("x", sprintf("%02d", 1:10))){
  d = get(i)
  d[1] <-  as.data.frame(apply(d[1],1,split.names))
  assign(i,d)
}


#combining all the gene.file
gene.file <- NULL
for(i in paste0("x", sprintf("%02d", 1:10))){
  d = get(i)
  gene.file = rbind(gene.file,d)
  assign(i,d)
}







