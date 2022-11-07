# Making genomic ranges for gwas and pre processing steps
preGWAS_process <- function(query.snp.gwas, n, organism){
  
  for(i in 1:n){ #sprintf("%02d", 1:n)
    x = query.snp.gwas
    d = get(paste0(x,i))
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
    #print(genename01)
    f <- get(paste0("zstat",i))
    f[,1]<- get(paste0("genename",i))
    f <- as.data.frame(t(f)) 
    colnames(f) <- f[1,]
    f <- f[-1,]
    f <- f %>% mutate_if(is.character,as.numeric, na.rm = T)
    assign(paste0("zstat",i),f)
    assign(i,f)
    g <- get(paste0("Marker",i))
    g[,1]<- get(paste0("genename",i))
    g <- as.data.frame(t(g)) 
    colnames(g) <- g[1,]
    g <- g[-1,]
    assign(paste0("Marker",i),g)
    assign(i,g)
    h <- get(paste0("pvalue",i))
    h[,1]<- get(paste0("genename",i))
    h <- as.data.frame(t(h)) 
    colnames(h) <- h[1,]
    h <- h[-1,]
    h <- h %>% mutate_if(is.character,as.numeric, na.rm = T)
    assign(paste0("pvalue",i),h)
    assign(i,h)
  }
}

preGWAS_process(query.snp.gwas,10)



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
  #print(genename01)
  f <- get(paste0("zstat",i))
  f[,1]<- get(paste0("genename",i))
  f <- as.data.frame(t(f)) 
  colnames(f) <- f[1,]
  f <- f[-1,]
  f <- f %>% mutate_if(is.character,as.numeric, na.rm = T)
  assign(paste0("zstat",i),f)
  #print(f[1:3,1:3])
  assign(i,f)
  g <- get(paste0("Marker",i))
  g[,1]<- get(paste0("genename",i))
  g <- as.data.frame(t(g)) 
  colnames(g) <- g[1,]
  g <- g[-1,]
  assign(paste0("Marker",i),g)
  #print(g[1:3,1:3])
  assign(i,g)
  h <- get(paste0("pvalue",i))
  h[,1]<- get(paste0("genename",i))
  h <- as.data.frame(t(h)) 
  colnames(h) <- h[1,]
  h <- h[-1,]
  h <- h %>% mutate_if(is.character,as.numeric, na.rm = T)
  assign(paste0("pvalue",i),h)
  #print(h[1:3,1:3])
  assign(i,h)
}

write.table(get(paste0("common",i)), paste0("common",i,".csv"), row.names = FALSE )
assign(paste0("common",i), vroom(paste0("common",i,".csv")))
system(paste0("rm common",i,".csv"))



# Subsetting chromosomes
setwd("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/All.phosphorus_LM")
folder <- getwd()
file_list <- list.files(path = folder, pattern = "*.txt")

# Reading all the files
for(i in 1:length(file_list)){
  assign(gsub(".txt","",file_list[i]), vroom(file_list[i]))
}


#Saving all the files
for(i in 1:length(file_list)){ #
  d = get(file_list[i])
  for(j in 1:10){
    x = paste0(file_list[i],"_",j)
    dplyr::filter(d, chr == j) %>% assign(x, ., inherits = TRUE)
    write.table(get(x),paste0(file_list[i],"_",j,".txt"), quote= FALSE, row.names = FALSE )
  } 
  assign(file_list[i], d)
}

# Remove extra .txt in the middle using shell
# > for file in *.txt
# > do
# > mv "$file" "${file/.txt/}"
# > done

system("sudo apt install rename")

#need to change the assign values and stuff mf
preGWAS_process(TP2.txt_,10)


