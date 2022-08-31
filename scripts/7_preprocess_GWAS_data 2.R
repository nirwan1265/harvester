###Data prep
##Filtering out the gene, fstat, Marker and pvalue columns
##Making a new dummy table
for(i in sprintf("%02d", 1:10)){
  assign(paste0("gwas",i), get(paste0("common",i)))
}

# #Saving RDS
# j <- 1
# for(i in paste0("common", sprintf("%02d", 1:10))){
#   d = get(i)
#   saveRDS(d, paste0("common",sprintf("%02d" , j),".RDS"))
#   assign(i,d)
#   j <- j+1
# }



##Filter table having only gene
for(i in paste0("gwas", sprintf("%02d", 1:10))){
  d=get(i)
  d <- d[which(d$first.X.Region == "gene"), ]
  assign(i,d)
}

##Filtering out the required columns for analysis
for(i in sprintf("%02d", 1:10)){
  assign(paste0("gwas",i), data.frame(get(paste0("gwas",i)))[,c(7,15,16,17)])
}

##Renaming columns
for(i in paste0("gwas", sprintf("%02d", 1:10))){
  d=get(i)
  colnames(d) = c("Gene","fstat","Marker","pvalue")
  assign(i,d)
}

##Sorting a/c gene name
#Not required, but might be useful for other database or when combining all results
# for(i in paste0("gwas", sprintf("%02d", 1:10))){
#   d=get(i)
#   d <- d[gtools::mixedorder(d$Gene), ]
#   assign(i,d)
# }

##Table with Z-stat values
for(i in paste0("gwas",sprintf("%02d",1:10))){
  d=get(i)
  assign(paste0(i,".zstat"), dcast(setDT(d), Gene~rowid(Gene, prefix = "zstat"), value.var = "zstat"))
  assign(i,d)
}
#----> continue
##Adding gene names
#Table with all the gene names from all the chromosome
#Can use this table for all the other tables with Marker and pvalue data
for(i in paste0("gwas", sprintf("%02d", 1:10))){
  d=get(i)
  assign(paste0(i,".gene.names"), get(paste0(i,".zstat"))[,1])
  assign(i,d)
}

for(i in paste0("gwas", sprintf("%02d", 1:10),".gene.names")){
  d=get(i)
  d <-  as.data.frame(apply(d,1,split.names))
  assign(i,d)
}

#----> continue
#Adding gene names
j <- 1
for(i in paste0("gwas",sprintf("%02d",1:10),".zstat")){
  d <- get(i)
  d[,1] <- get(paste0("gwas",sprintf("%02d",j),".gene.names"))[,1]
  assign(i,d)
  j = j + 1
}
#Data ordering
for(i in paste0("gwas",sprintf("%02d", 1:10),".zstat")){
  d <- get(i)
  d <- as.data.frame(t(d))
  x <- d[1,]
  d <- d[-1,]
  d <- apply(d,2,sort)
  d <- as.data.frame(stri_list2matrix(d, byrow = FALSE))
  colnames(d) <- x
  assign(i,d)
}


##Table with SNP markers
for(i in paste0("gwas",sprintf("%02d",1:10))){
  d=get(i)
  assign(paste0(i,".Marker"), dcast(setDT(d), Gene~rowid(Gene, prefix = "Marker"), value.var = "Marker"))
  assign(i,d)
}
#Adding gene names
j <- 1
for(i in paste0("gwas",sprintf("%02d",1:10),".Marker")){
  d <- get(i)
  d[,1] <- get(paste0("gwas",sprintf("%02d",j),".gene.names"))[,1]
  assign(i,d)
  j = j + 1
}
#Data ordering
for(i in paste0("gwas",sprintf("%02d", 1:10),".Marker")){
  d <- get(i)
  d <- as.data.frame(t(d))
  x <- d[1,]
  d <- d[-1,]
  d <- apply(d,2,sort)
  d <- as.data.frame(stri_list2matrix(d, byrow = FALSE))
  colnames(d) <- x
  assign(i,d)
}

##Table with pvalue values
for(i in paste0("gwas",sprintf("%02d",1:10))){
  d=get(i)
  assign(paste0(i,".pvalue"), dcast(setDT(d), Gene~rowid(Gene, prefix = "pvalue"), value.var = "pvalue"))
  assign(i,d)
}
#Adding gene names
j <- 1
for(i in paste0("gwas",sprintf("%02d",1:10),".pvalue")){
  d <- get(i)
  d[,1] <- get(paste0("gwas",sprintf("%02d",j),".gene.names"))[,1]
  assign(i,d)
  j = j + 1
}
#Data ordering
for(i in paste0("gwas",sprintf("%02d", 1:10),".pvalue")){
  d <- get(i)
  d <- as.data.frame(t(d))
  x <- d[1,]
  d <- d[-1,]
  d <- apply(d,2,sort)
  d <- as.data.frame(stri_list2matrix(d, byrow = FALSE))
  colnames(d) <- x
  assign(i,d)
}


#saving objects
# j <- 1
# for(i in paste0("gwas", sprintf("%02d", 1:10),".fstat")){
#   d = get(i)
#   saveRDS(d, paste0("gwas",sprintf("%02d" , j),".fstat.RDS"))
#   assign(i,d)
#   j <- j+1
# }
# j <- 1
# for(i in paste0("gwas", sprintf("%02d", 1:10),".gene.names")){
#   d = get(i)
#   saveRDS(d, paste0("gwas",sprintf("%02d" , j),".gene.names.RDS"))
#   assign(i,d)
#   j <- j+1
# }
# j <- 1
# for(i in paste0("gwas", sprintf("%02d", 1:10),".Marker")){
#   d = get(i)
#   saveRDS(d, paste0("gwas",sprintf("%02d" , j),".Marker.RDS"))
#   assign(i,d)
#   j <- j+1
# }
# j <- 1
# for(i in paste0("gwas", sprintf("%02d", 1:10),".pvalue")){
#   d = get(i)
#   saveRDS(d, paste0("gwas",sprintf("%02d" , j),".pvalue.RDS"))
#   assign(i,d)
#   j <- j+1
# }
# 


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

