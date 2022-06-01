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

#phenotype:



