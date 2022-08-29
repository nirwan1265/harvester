#Global Functions:
## Database Annotation
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Maize.annotation")
db <- as.data.frame(read.table(file ="Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3", sep = "\t", header = FALSE))
dbannot <- function(n,db,output){
  colnames(db) <- c("Chromosome","Database","Region","Start","End","NA","Strand","NA2","Gene")
  if(db[1,1] == "chr1"){
    for(i in 1:n){
      assign(paste0("db.chr.",i), db[which(db$Chromosome == paste0("chr",i)), ])
    }
  }
  if(db[1,1] == "Chr1"){
    for(i in 1:n){
      assign(paste0("db.chr.",i), db[which(db$Chromosome == paste0("Chr",i)), ])
    }
  }
  
  return(as.data.frame(db.chr.1))
}


##Gene Name filtering
split.names <- function(x,split){
  split.genename <- unlist(strsplit(x, split = ';', fixed = TRUE))[1]
  split.genename2 <- unlist(strsplit(split.genename, split = ":", fixed = TRUE))[2]
  return(split.genename2)
}

## Converting pvalues to Z-scores
zval <- function(x, output){
  pvalue <- unlist(as.numeric(x[4]))
  o <- p.to.Z(pvalue)
  return(o)
}

# ACAT function
##For combining known number of pvalues, replace 2 in 1:2 by the number of pvalues
acat <- function(x, output){
  pvalue <- as.numeric(x[1:2])
  pvalue <- unlist(pvalue)
  o <- ACAT(Pvals = pvalue)
  return(o)
}
## Usage
### dataframe <- cbind(pvalue.dataframe, (apply(pvalue.dataframe[,c(1:n),1,acat])))

## For combining n number of pvalues
acat <- function(x, n, output){
  pvalue <- as.numeric(x[1:n])
  pvalue <- unlist(pvalue)
  o <- ACAT(Pvals = pvalue)
  return(o)
}

## Usage
### as.data.frame(apply(pvalue.dataframe,1,acat, n = <number of pvalues in the columns>))

##Pvalues Combinations Test
pvalue.combine <- function(gwas.zstat, gwas.marker, gwas.pvalue, geno, tab.pc,combined.test.statistics){
  x <- as.data.frame(matrix(0, nrow = 1, ncol = 1))
  y <- vector()
  z <- vector()
  combined.test.statistics <- as.data.frame(matrix(NA, nrow = 1, ncol = 1))
  ref_genotype <- as.data.frame(matrix(NA, nrow = 1, ncol = 1))
  ref_genotype_skat <- as.data.frame(matrix(NA, nrow = 1, ncol = 1))
  
  for (i in 1:ncol(gwas.zstat)){ 
    for(j in 1:sum(!is.na(gwas.zstat[,i]))){ 
      x[1,j] <- gwas.zstat[j,i]
      
      y[j] <- as.vector(as.character(gwas.marker[j,i]))
      
      z[j] <- gwas.pvalue[j,i]
      z <- as.double(z[!is.na(z)])
    }
    x <- as.matrix(as.double(x))
    if(nrow(x) > 2000){
      x <- x[1:2000,]
      y <- y[1:2000]
      z <- z[1:2000]
      
      #GBJ, minP, GHC, OMNI
      ref_genotype <- as.data.frame(geno[,colnames(geno) %in% y])
      # ref_genotype <- data.frame(lapply(ref_genotype, function(x){
      #   gsub("-",0,x)
      # }))
      # ref_genotype <- data.frame(apply(ref_genotype, 2, function(x) as.numeric(as.character(x))))
      cor_mat <- estimate_ss_cor(ref_pcs=tab.pc, ref_genotypes=ref_genotype, link_function='linear')
      gbj.test <- GBJ(test_stats = x, cor_mat=cor_mat)
      ghc.test <- GHC(test_stats = x, cor_mat=cor_mat)
      minP.test <- minP(test_stats = x, cor_mat=cor_mat)
      OMNI.test <- OMNI_ss(test_stats = x, cor_mat=cor_mat, num_boots = 100)
      combined.test.statistics[i,1] <- gbj.test$GBJ_pvalue
      combined.test.statistics[i,2] <- ghc.test$GHC_pvalue
      combined.test.statistics[i,3] <- minP.test$minP_pvalue
      
      #SKAT
      ref_genotype_skat <- as.matrix(ref_genotype)
      obj01 <- as.list(ref_genotype_skat,pheno)
      obj01 <- SKAT_Null_Model(pheno ~ 1, out_type="C", data=obj01)
      combined.test.statistics[i,4] <- SKAT(ref_genotype_skat,obj01)$p.value
      
      #OMNI
      combined.test.statistics[i,5] <- OMNI.test$OMNI_pvalue
      
      
      x <- as.data.frame(matrix(0, nrow = 1, ncol = 1))
      y <- vector()
      
    } else if(nrow(x) >= 2 & nrow(x) < 2000){
      ref_genotype <- as.data.frame(geno[,colnames(geno) %in% y])
      cor_mat <- estimate_ss_cor(ref_pcs=tab.pc, ref_genotypes=ref_genotype, link_function='linear')
      
      #GBJ, minP, GHC, OMNI
      gbj.test <- GBJ(test_stats = x, cor_mat=cor_mat)
      minP.test <- minP(test_stats = x, cor_mat=cor_mat)
      ghc.test <- GHC(test_stats = x, cor_mat=cor_mat)
      OMNI.test <- OMNI_ss(test_stats = x, cor_mat=cor_mat, num_boots = 100)
      combined.test.statistics[i,1] <- gbj.test$GBJ_pvalue
      combined.test.statistics[i,2] <- ghc.test$GHC_pvalue
      combined.test.statistics[i,3] <- minP.test$minP_pvalue
      
      #SKAT
      ref_genotype_skat <- as.matrix(ref_genotype)
      obj01 <- as.list(ref_genotype_skat,pheno)
      obj01 <- SKAT_Null_Model(pheno ~ 1, out_type="C", data=obj01)
      combined.test.statistics[i,4] <- SKAT(ref_genotype_skat,obj01)$p.value
      
      #OMNI
      combined.test.statistics[i,5] <- OMNI.test$OMNI_pvalue
      
      x <- as.data.frame(matrix(0, nrow = 1, ncol = 1))
      y <- vector()
      ref_genotype <- as.data.frame(matrix(NA, nrow = 1, ncol = 1))
      ref_genotype_skat <- as.data.frame(matrix(NA, nrow = 1, ncol = 1))
      
    } else if(nrow(x) == 1){
      combined.test.statistics[i,1] <- as.double(gwas.pvalue[1,i])
      combined.test.statistics[i,2] <- as.double(gwas.pvalue[1,i])
      combined.test.statistics[i,3] <- as.double(gwas.pvalue[1,i])
      combined.test.statistics[i,4] <- as.double(gwas.pvalue[1,i])
      combined.test.statistics[i,5] <- as.double(gwas.pvalue[1,i])
      
      x <- as.data.frame(matrix(0, nrow = 1, ncol = 1))
      y <- vector()
      ref_genotype <- as.data.frame(matrix(NA, nrow = 1, ncol = 1))
      ref_genotype_skat <- as.data.frame(matrix(NA, nrow = 1, ncol = 1))
    }
    #combined.test.statistics[i,6] <- apply(combined.test.statistics[i,1:4],1,acat)
    x <- as.data.frame(matrix(0, nrow = 1, ncol = 1))
    y <- vector()
    z <- vector()
    ref_genotype <- as.data.frame(matrix(NA, nrow = 1, ncol = 1))
    ref_genotype_SKAT <- as.data.frame(matrix(NA, nrow = 1, ncol = 1))
  }
  #colnames(combined.test.statistics) <- c("GBJ","GHC","minP","SKAT","OMNI(1-4)","CCT(1-4)")
  colnames(combined.test.statistics) <- c("GBJ","GHC","minP","SKAT","OMNI(1-4)")
  return(combined.test.statistics)
}

