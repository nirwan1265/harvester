################################################################################
#SNP FILE
################################################################################

#Collecting all the genotype SNP names
for(i in sprintf("%02d", 1:10)){
  assign(paste0("x",i), as.data.frame(colnames(get(paste0("geno", i)))))
}

#Chromosome Name
j <- 1 
for(i in paste0("x", sprintf("%02d", 1:10))){
  d = get(i)
  d$chr = paste0("chr",j)
  j = j + 1
  assign(i,d)
}

#Combining all the SNPs
x <- NULL
for(i in paste0("x", sprintf("%02d", 1:10))){
  d = get(i)
  x = rbind(x,d[1])
  assign(i,d)
}

#Combining all the chr
y <- NULL
for(i in paste0("x", sprintf("%02d", 1:10))){
  d = get(i)
  y = rbind(y,d[2])
  assign(i,d)
}


#Combining SNPs and Chr
snp.file <- cbind(x,y)

#Function for extracting the base pair
base.pair <- function(x,split){
  bp <- unlist(strsplit(x, split = '_', fixed = TRUE))[2]
  return(bp)
}

for(i in paste0("x", sprintf("%02d", 1:10))){
  d=get(i)
  d <-  as.data.frame(apply(d,1,base.pair))
  assign(i,d)
}


#Combining all the base pairs
z <- NULL
for(i in paste0("x", sprintf("%02d", 1:10))){
  d = get(i)
  z = rbind(z,d[1])
  assign(i,d)
}

#Combining SNPS, Chr, and BP
snp.file <- cbind(gene.file, z)

colnames(snp.file) <- c("SNP ID", "chromosome", "base pair position")

#query.snp has it all


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








