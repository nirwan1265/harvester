#Raw GWAS result
x <- gwas10.pvalue
x[] <- lapply(x, function(x) as.numeric(replace(x, is.na(x), 1)))

raw.gwas <- as.data.frame(apply(x, 2, min))
colnames(raw.gwas) <- "pvalue"

raw.gwas <- as.data.frame(raw.gwas[order(raw.gwas$pvalue), ])



#OMNI GWAS result

