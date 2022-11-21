GWAS_process <- function(path, filename, n, organism){
  if(organism == "Sorghum bicolor"){
    file_list <- list.files(path = path, pattern = filename)
    for(i in 1:length(file_list)){
      #assign(gsub(".txt","",file_list[i]), vroom(file_list[i]))
      assign(paste0(file_list[i],"sorghum"), vroom(file_list[i]))
      #write.table(file_list[i], paste0(file_list[i],".csv"), row.names = FALSE )
      #assign(paste0(file_list[i]), read.csv(paste0(file_list[i],".csv")))
      #system(rm(paste0(file_list[i],".csv")))
      #preGWAS_process(file_list[i])
      #print(get(file_list[i]))
      #print(tot1.txt)
      #preGWAS_process(get(file_list[i]))
      #assign(file_list[i])
      preGWAS_process(paste0(file_list[i],"sorghum"))
      
    }
    #return(file_list)
    
    #assign(file_list[i])
  }
  # elseif(organism == "Zea mays"){
  # print("hello")  
  # }
}




preGWAS_process <- function(query.snp.gwas){
  a <- 1
  for(i in query.snp.gwas){
    d = as.data.frame(get(i))
    print("d is",d)
    assign(i,d)
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
    assign(i,d,e,f,g,h)
    
    # write.table(get(paste0("common",i)), paste0("common",i,".csv"), row.names = FALSE )
    # assign(paste0("common",i), vroom(paste0("common",i,".csv")))
    # system(paste0("rm common",i,".csv"))
    
  }
}

setwd("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/All.phosphorus_LM")
path <- getwd()
filename <- "tot"
organism <- "Sorghum bicolor"
GWAS_process(path, filename, 10,  organism)
system("ls")




preprocess <- function(path, filename, n, organism){
  a <- 1
  if(organism == "Sorghum bicolor"){
    file_list <- list.files(path = path, pattern = filename)
    for(i in 1:length(file_list)){
      assign(file_list[i], vroom(file_list[i]))
      d = as.data.frame(get(file_list[i]))
      d$zstat = unlist(apply(d,1,zval))
      names(d) <- c("Marker","chr","Start_Position","pvalue","Zvalue")
      d <- d %>% mutate_at(c('chr','Start_Position','pvalue','Zvalue'),as.numeric)
      #print(d[1:5,1:5])
      assign(paste0("gr.q", i) , GRanges(seqnames = paste0("chr",i), ranges = IRanges(start = d[,"Start_Position"], width = 1, zstat = d[,"Zvalue"], Marker = d[,"Marker"],pvalue = d[,"pvalue"])))
      assign(paste0("common",i), as.data.frame(findOverlapPairs(get(paste0("gr.db",sprintf("%02d",a))), get(paste0("gr.q",sprintf("%02d",a))))))
      a = a + 1
      e = get(paste0("common",i))
      e = e[which(e$first.X.Region == "gene"), ]
      e = e[,c(7,16,17,18)]
      colnames(e) = c("Gene","Marker","pvalue","zstat")
      #print(e[1:4,1:4])
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
      #print(f[1:3,1:3])
      assign(paste0("zstat",i),f)
      #assign(i,f)
      g <- get(paste0("Marker",i))
      g[,1]<- get(paste0("genename",i))
      g <- as.data.frame(t(g))
      colnames(g) <- g[1,]
      g <- g[-1,]
      #print(g[1:3,1:3])
      assign(paste0("Marker",i),g)
      #assign(i,g)
      h <- get(paste0("pvalue",i))
      h[,1]<- get(paste0("genename",i))
      h <- as.data.frame(t(h))
      colnames(h) <- h[1,]
      h <- h[-1,]
      h <- h %>% mutate_if(is.character,as.numeric, na.rm = T)
      print(h[1:3,1:3])
      assign(paste0("pvalue",i),h)
      #assign(i,d,e,f,g,h)
    }
    
    #assign(d)
    #return(file_list[i])
  }
  else if(organism == "Zea mays"){
    print("hello")
  }
}


setwd("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/All.phosphorus_LM")
path <- getwd()
filename <- "tot"
organism <- "Sorghum bicolor"
precomb(path, filename, 10,  organism)

