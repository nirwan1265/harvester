###Annotation the SNPs
##Making GRanges for Database 
for(i in sprintf("%02d", 1:10)){
  assign(paste0("gr.db", i) , GRanges(seqnames = paste0("chr",i), ranges = IRanges(start = get(paste0("db.",i))[,"Start"], end = get(paste0("db.",i))[,"End"]), strand = get(paste0("db.",i))[,"Strand"], Region = get(paste0("db.",i))[,"Region"], Gene = get(paste0("db.",i))[,"Gene"]))
}

## for maize
for(i in sprintf("%02d", 1:10)){
  assign(paste0("gr.db", i) , GRanges(seqnames = paste0("chr",i), ranges = IRanges(start = get(paste0("ref_",i))[,"Start"], end = get(paste0("ref_",i))[,"End"]), strand = get(paste0("ref_",i))[,"Strand"], Region = get(paste0("ref_",i))[,"Region"], Gene = get(paste0("ref_",i))[,"Gene"]))
}
a <- 1
for(i in paste0("ref_", sprintf("%02d", 1:10))){
  saveRDS(i, paste0("gr.db", sprintf("%02d",a),".RDS"))
  a <- a+1
}



##Making GRanges for gwas Query
#Changing the position column to numeric and removing the first row
for(i in paste0("query.snp.gwas", sprintf("%02d", 1:10))){
  d = get(i)
  d = d[-1,]
  d$Pos = as.numeric(d$MLM_Stats.Pos)
  assign(i,d)
}

for(i in sprintf("%02d", 1:10)){
  assign(paste0("gr.q", i) , GRanges(seqnames = paste0("chr",i), ranges = IRanges(start = get(paste0("query.snp.gwas",i))[,"Pos"], width = 1, fstat = get(paste0("query.snp.gwas",i))[,"MLM_Stats.F"], Marker = get(paste0("query.snp.gwas",i))[,"MLM_Stats.Marker"],pvalue = get(paste0("query.snp.gwas",i))[,"MLM_Stats.p"])))
}


#Finding the Overlaps
for(i in sprintf("%02d", 1:10)){
  assign(paste0("common",i), as.data.frame(findOverlapPairs(get(paste0("gr.db",i)), get(paste0("gr.q",i)))))
}




