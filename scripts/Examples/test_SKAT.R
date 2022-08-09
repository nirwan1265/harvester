library(SKAT)
data("SKAT.example")
attributes(SKAT.example)

#Covariates
X <- SKAT.example$X
X

#Continuous phenotype
y.c <- SKAT.example$y.c

#binaru phenotype
y.b <- SKAT.example$y.b
y.b
#############################################################
# SKAT with default Beta(1,25) Weights
# - without covariates
Z<-SKAT.example$Z

# continuous trait
obj<-SKAT_Null_Model(y.c ~ 1, out_type="C", data=SKAT.example)
SKAT(Z, obj)$p.value




##################################################
# SKAT with default Beta(1,25) Weights
# - Optimal Test
SKAT(Z, obj, method="optimal.adj")$p.value
# you can get the same p-value by using method="SKATO"
SKAT(Z, obj, method="SKATO")$p.value



#############################################################
# SKAT with Beta(1,30) Weights
SKAT(Z, obj, weights.beta=c(1,30))$p.value



#############################################################
#Genotype
setwd("/Users/nirwan/Library/Mobile Documents/com~apple~CloudDocs/Data for sorghum/sorghum/Lasky.hapmap/90perc.taxa.filter.95perc.site.filter/numerical")
system("ls")
geno01 <- read.delim("chr01.numeric.txt")
geno01.trial <- geno01[,1:10]
geno01.trial <- as.matrix(geno01.trial)

#phenotype:
setwd("/Users/nirwan/Library/Mobile Documents/com~apple~CloudDocs/Data for sorghum/sorghum/Lasky.hapmap/filtered.phenotype")
pheno01 <- read.table("phenotype.chr01.txt", header = FALSE)
row.names(pheno01) <- pheno01[,1]
pheno01 <- pheno01[,-1]
pheno01 <- as.matrix(pheno01)


#Creating object
obj01 <- as.list(geno01.trial,pheno01)
obj01<-SKAT_Null_Model(pheno01 ~ 1, out_type="C", data=obj01)

#pvalue
SKAT(geno01.trial, obj01)$p.value



typeof(Z)
class(Z)


####
setwd("/Users/nirwan/Library/Mobile Documents/com~apple~CloudDocs/Data for sorghum/sorghum/Lasky.hapmap/90perc.taxa.filter.95perc.site.filter/phenotype")
pheno01 <- read.csv("phenotype.chr01.csv", header = FALSE)
pheno01 <- as.matrix(pheno01)

x <- as.data.frame(matrix(0, nrow = 1, ncol = 1))
y <- vector()
z <- vector()
combined.test.statistics <- as.data.frame(matrix(NA, nrow = 1, ncol = 1))
ref_genotype <- as.data.frame(matrix(NA, nrow = 1, ncol = 1))

for (i in 10){ #ncol(gwas1.Test.Stat)
  for(j in 1:sum(!is.na(gwas01.fstat[,i]))){ 
    x[1,j] <- gwas01.fstat[j,i]
    #x <- as.double(x[!is.na(x)])
    
    y[j] <- as.vector(as.character(gwas01.Marker[j,i]))
    #y <- y[!is.na(y)]
    
    z[j] <- gwas01.pvalue[j,i]
    z <- as.double(z[!is.na(z)])
  }
  x <- as.matrix(as.double(x))
  if(nrow(x) > 2000){
    x <- x[1:2000,]
    y <- y[1:2000]
    z <- z[1:2000]
    ref_genotype <- as.data.frame(geno1[,colnames(geno01) %in% y])
    ref_genotype <- data.frame(lapply(ref_genotype, function(x){
      gsub("-",0,x)
    }))
    ref_genotype <- data.frame(apply(ref_genotype, 2, function(x) as.numeric(as.character(x))))
    cor_mat <- estimate_ss_cor(ref_pcs=tab.pc01, ref_genotypes=ref_genotype, link_function='linear')
    bj.test <- BJ(test_stats = x, cor_mat=cor_mat)
    gbj.test <- GBJ(test_stats = x, cor_mat=cor_mat)
    minP.test <- minP(test_stats = x, cor_mat=cor_mat)
    hc.test <- HC(test_stats = x, cor_mat=cor_mat)
    ghc.test <- GHC(test_stats = x, cor_mat=cor_mat)
    combined.test.statistics[i,1] <- bj.test$BJ_pvalue
    combined.test.statistics[i,2] <- gbj.test$GBJ_pvalue
    combined.test.statistics[i,3] <- hc.test$HC_pvalue
    combined.test.statistics[i,4] <- ghc.test$GHC_pvalue
    combined.test.statistics[i,5] <- minP.test$minP_pvalue
    x <- as.data.frame(matrix(0, nrow = 1, ncol = 1))
    y <- vector()
  } else if(nrow(x) >= 2 & nrow(x) < 2000){
    ref_genotype <- as.data.frame(geno01[,colnames(geno01) %in% y])
    ref_genotype <- data.frame(lapply(ref_genotype, function(x){
      gsub("-",0,x)
    }))
    ref_genotype <- data.frame(lapply(ref_genotype, function(x){
      gsub("0.5",9,x)
    }))
    ref_genotype <- data.frame(apply(ref_genotype, 2, function(x) as.numeric(as.character(x))))
    ref_genotype <- as.matrix(ref_genotype)
    
    obj01 <- as.list(ref_genotype,pheno01)
    obj01<-SKAT_Null_Model(pheno01 ~ 1, out_type="C", data=obj01)
    
    combined.test.statistics[i,1] <- SKAT(ref_genotype,obj01)$p.value
    
    x <- as.data.frame(matrix(0, nrow = 1, ncol = 1))
    y <- vector()
  } else if(nrow(x) == 1){
    combined.test.statistics[i,1] <- as.double(gwas01.pvalue[1,i])
  }
  ref_genotype <- as.data.frame(matrix(NA, nrow = 1, ncol = 1))
  x <- as.data.frame(matrix(0, nrow = 1, ncol = 1))
  y <- vector()
  z <- vector()
}

x
y
ref_genotype

system("ls")
setwd("/Users/nirwan/Library/Mobile Documents/com~apple~CloudDocs/Data for sorghum/sorghum/Lasky.hapmap/")
