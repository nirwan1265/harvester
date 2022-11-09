preGWAS_process <- function(path, filename, n, organism){
  if(organism == "Sorghum bicolor"){
    file_list <- list.files(path = path, pattern = filename)
    for(i in 1:length(file_list)){
      #assign(gsub(".txt","",file_list[i]), vroom(file_list[i]))
      assign(file_list[i], vroom(file_list[i]))
     #GWAS_process(file_list[i])
    }
    #return(file_list[i])
    GWAS_process(file_list[i])
    #assign(file_list[i],i)
  }
  # elseif(organism == "Zea mays"){
  # print("hello")  
  # }
}

GWAS_process <- function(query.snp.gwas){
  a <- 1
  for(i in query.snp.gwas){
    d = get(i)
    print(d)
    assign(i,d)
  }
}

preGWAS_process(path, filename, 10,  organism)


GWAS_process <- function(query.snp.gwas){
  a <- 1
  for(i in query.snp.gwas){
    d = get(i)
    print(d)
    d$zstat = unlist(apply(d,1,zval))
    names(d) <- c("Marker","chr","Start_Position","pvalue","Zvalue")
    d <- d %>% mutate_at(c('chr','Start_Position','pvalue','Zvalue'),as.numeric)
    print(d[1:5,1:5])
    assign(paste0("gr.q", i) , GRanges(seqnames = paste0("chr",i), ranges = IRanges(start = d[,"Start_Position"], width = 1, zstat = d[,"Zvalue"], Marker = d[,"Marker"],pvalue = d[,"pvalue"])))
    assign(paste0("common",i), as.data.frame(findOverlapPairs(get(paste0("gr.db",sprintf("%02d",a))), get(paste0("gr.q",sprintf("%02d",a))))))
    a = a + 1
    e = get(paste0("common",i))
    e = e[which(e$first.X.Region == "gene"), ]
    e = e[,c(7,16,17,18)]
    colnames(e) = c("Gene","zstat","Marker","pvalue")
    print(e[,1:2])
    assign(paste0("filter_common", i), e)
    assign(paste0("zstat",i), dcast(setDT(e), Gene~rowid(Gene, prefix = "zstat"), value.var = "zstat"))
    assign(paste0("Marker",i), dcast(setDT(e), Gene~rowid(Gene, prefix = "Marker"), value.var = "Marker"))
    assign(paste0("pvalue",i), dcast(setDT(e), Gene~rowid(Gene, prefix = "pvalue"), value.var = "pvalue"))
    assign(paste0("genename",i),apply(get(paste0("Marker",i)),1,split.names))
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
    assign(i,d,e)
  }
}

setwd("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/All.phosphorus_LM")
path <- getwd()
filename <- "tot"
organism <- "Sorghum bicolor"
preGWAS_process(path, filename, 10,  organism)


