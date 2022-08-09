library(devtools)
devtools::install_github("yaowuliu/ACAT")
library(ACAT)


#Reading files
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/pvalues.combination")
pvalue.chr01 <- read.csv("combined.omni.magma.csv")
rownames(pvalue.chr01) <- pvalue.chr01[,1]
pvalue.chr01 <- pvalue.chr01[,-1]

class(pvalue.chr01)

#Aggregated Cauchy Association Test
pvalue <- as.numeric(pvalue.chr01[1,])
pvalue
typeof(pvalue)
class(pvalue)
ACAT(Pvals = pvalue)
ACAT(matrix(runif(1000), ncol = 10))

#pvalues equal to 1
#A: When a test statistic follows a continuous ditribution, the corresponding p-value should follow a uniform distribution between 0 and 1 under the null. Hence, for continuous distributions, one should never get a p-value that is exactly 1. In practice, we may have p-values being 1 becasue the calibrated p-value is an approximation of the "true" p-value. For example, if simulation-based methods (e.g., permutation) is used, one could have p-values equal to 1 since the number of simulations is always finite.
#Therefore, if there are p-values equal to 1, we need to first find out the reason and then try to "correct" the calibrated p-values. For example, if it is due to permutation, we can replace 1 by 1-1/N, where N is the number of permutations. Another simple way is to replace 1 by 1-1/d, where d is the number of p-values combined by ACAT.

#ACAT function
acat <- function(x, output){
  pvalue <- unlist(as.numeric(x[1]))
  o <- ACAT(Pvals = pvalue)
  return(o)
}

y <- chr01[,-1]
x <-as.data.frame(apply(y,1,acat))
pvalue <- chr01[]
y

chr01 <- 
