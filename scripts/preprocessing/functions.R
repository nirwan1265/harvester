# Making genomic ranges for gwas results
grangeGWAS <- function(query.snp.gwas, n){
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
    print(assign(paste0("common",i), as.data.frame(findOverlapPairs(get(paste0("gr.db",i)), get(paste0("gr.q",i))))))
    # x <- get(paste0("common",i))
    # y <- rbind(x,x)
    # return(y)
  }
}

x <- grangeGWAS(query.snp.gwas,10)


grangeGWAS <- function(query.snp.gwas, n, organism){
  for(i in sprintf("%02d", 1:10)){
    d = get(paste0("query.snp.gwas", i))
    d = d[-1,]
    #d <-  as_tibble(d)
    names(d) <- c("Marker","chr","Start_Position","Zvalue","pvalue","End_Position")
    d <- d %>% mutate_at(c('chr','Start_Position','Zvalue','pvalue','End_Position'),as.integer)
    assign(paste0("query.gwas", i), d)
    assign(paste0("gr.q", i) , GRanges(seqnames = paste0("chr",i), ranges = IRanges(start = get(paste0("query.gwas",i))[,"Start_Position"], width = 1, zstat = get(paste0("query.gwas",i))[,"Zvalue"], Marker = get(paste0("query.gwas",i))[,"Marker"],pvalue = get(paste0("query.gwas",i))[,"pvalue"])))
    assign(paste0("common",i), as.data.frame(findOverlapPairs(get(paste0("gr.db",i)), get(paste0("gr.q",i)))))
    write.csv(get(paste0("common",i)), paste0("common",i,".csv"), row.names = FALSE)
    assign(paste0("common",i), vroom(paste0("common",i,".csv")))
    system(paste0("rm common",i,".csv"))
    assign(i,d)
  }
}

grangeGWAS(query.snp.gwas,10)



for(i in sprintf("%02d", 1:10)){
  d = get(paste0("query.snp.gwas", i))
  d = d[-1,]
  names(d) <- c("Marker","chr","Start_Position","Zvalue","pvalue")
  d <- d %>% mutate_at(c('chr','Start_Position','Zvalue','pvalue'),as.numeric)
  assign(paste0("query.gwas", i), d)
  assign(paste0("gr.q", i) , GRanges(seqnames = paste0("chr",i), ranges = IRanges(start = get(paste0("query.gwas",i))[,"Start_Position"], width = 1, zstat = get(paste0("query.gwas",i))[,"Zvalue"], Marker = get(paste0("query.gwas",i))[,"Marker"],pvalue = get(paste0("query.gwas",i))[,"pvalue"])))
  assign(paste0("common",i), as.data.frame(findOverlapPairs(get(paste0("gr.db",i)), get(paste0("gr.q",i)))))
  assign(i,d)
  e = get(paste0("common",i))
  e <- e[which(e$first.X.Region == "gene"), ]
  e = e[,c(7,15,16,17)]
  colnames(e) = c("Gene","zstat","Marker","pvalue")
  e$zstat = unlist(apply(e,1,zval))
  assign(paste0("filter_common", i), e)
  assign(paste0("zstat",i), dcast(setDT(e), Gene~rowid(Gene, prefix = "zstat"), value.var = "zstat"))
  assign(paste0("Marker",i), dcast(setDT(e), Gene~rowid(Gene, prefix = "Marker"), value.var = "Marker"))
  assign(paste0("pvalue",i), dcast(setDT(e), Gene~rowid(Gene, prefix = "pvalue"), value.var = "pvalue"))
  assign(paste0("genename",i),apply(get(paste0("Marker",i)),1,split.names))
  assign(i,e)
  print(genename01)
  f =  get(paste0("zstat",i))
  f[,1]<- get(paste0("genename",i))
  print(f)
  assign(i,f)
}

##Gene Name filtering
split.names <- function(x,split){
  split.genename <- unlist(strsplit(x, split = ';', fixed = TRUE))[1]
  split.genename2 <- unlist(strsplit(split.genename, split = ":", fixed = TRUE))[2]
  return(split.genename2)
}

Marker01<- apply(Marker01,1,split.names)
gene01

Marker01 <- replace(Marker01$Gene,gene01)

common01[1]
system("ls")
system("pwd")
rm(common01)

write.table(get(paste0("common",i)), paste0("common",i,".csv"), row.names = FALSE )
assign(paste0("common",i), vroom(paste0("common",i,".csv")))
system(paste0("rm common",i,".csv"))