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
# j = 1
# for(i in paste0("pvalue.combine",sprintf("%02d", 1:10))){
#   d = get(i)
#   a <- get(paste0("gwas",sprintf("%02d", j), ".gene.names"))
#   a <- unlist(a)
#   row.names(d) <- a
#   j = j + 1
#   assign(i,d)
# }


#Column names
# for(i in paste0("pvalue.combine",sprintf("%02d", 1:10))){
#   d = get(i)
#   colnames(d) <- c("GBJ","GHC","minP","SKAT","OMNI_ItoIV","CCT_ItoIV")
#   assign(i,d)
# }


#Reading the files
for(i in sprintf("%02d", 1:10)){
  assign(paste0("pvalue.combine",i), read.csv(paste0("pvalue.combine",i,".csv"), row.names = 1))
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


# USING ACAT
j <- 1
for(i in paste0("pvalue.combine", sprintf("%02d", 1:10))){
  d = get(i)
  assign(paste0("cauchy",sprintf("%02d",j)),as.data.frame(apply(d,1,acat, n = 4))) # n = number of rows
  j <- j + 1
  assign(i,d)
}

# Adding ACAT to Pvalue
pvalue.combine01$CCT_ItoIV <- unlist(cauchy01)
pvalue.combine02$CCT_ItoIV <- unlist(cauchy02)
pvalue.combine03$CCT_ItoIV <- unlist(cauchy03)
pvalue.combine04$CCT_ItoIV <- unlist(cauchy04)
pvalue.combine05$CCT_ItoIV <- unlist(cauchy05)
pvalue.combine06$CCT_ItoIV <- unlist(cauchy06)
pvalue.combine07$CCT_ItoIV <- unlist(cauchy07)
pvalue.combine08$CCT_ItoIV <- unlist(cauchy08)
pvalue.combine09$CCT_ItoIV <- unlist(cauchy09)
pvalue.combine10$CCT_ItoIV <- unlist(cauchy10)



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


#SAVING the files
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/pvalues.combination")
j = 1
for(i in paste0("pvalue.combine",sprintf("%02d", 1:10))){
  d = get(i)
  write.csv(d, paste0("pvalue.combine",sprintf("%02d", j),".csv"))
  j = j + 1
  assign(i,d)
}

