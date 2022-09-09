Z<-SKAT.example$Z
Z
obj<-SKAT_Null_Model(y.c ~ X, out_type="C", data=SKAT.example)
x <- SKAT.example

weights<-Get_Logistic_Weights(Z, par1=0.07, par2=150)
SKAT(Z, obj, kernel = "linear.weighted", weights=weights)$p.value

?SKAT_NULL_emmaX()


data=SKAT.fam.example
dim(SKAT.fam.example$K)

obj<-SKAT_NULL_emmaX(y ~ X, K=K, data=SKAT.fam.example)

K = SKAT.fam.example$K
Z = SKAT.fam.example$Z
y <- SKAT.fam.example$y
X <- SKAT.fam.example$X
typeof(y)
typeof(pheno)
class(y)
class(pheno)
typeof(Z)
typeof(ref_genotype_skat)
class(Z)
class(ref_genotype_skat)

nrow(K)
nrow(Z)
nrow(y)
nrow(X)
typeof(K)
typeof(k)
class(K)
class(k)
typeof(K[1,1])
typeof(k[1,1])

nrow(pheno)
nrow(ref_genotype)
nrow(cor_mat)
cor_mat


k <- kinship(
  as.matrix(ref_genotype),
  method = c("IBS"),
  denominator = NULL
)

typeof(ref_genotype)
class()

install.packages("statgenGWAS")
library(statgenGWAS)
