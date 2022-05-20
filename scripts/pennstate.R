#Raw GWAS result
setwd("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Data for sorghum/sorghum/Results")
x <- gwas01.pvalue
x[] <- lapply(x, function(x) as.numeric(replace(x, is.na(x), 1)))

raw.gwas <- as.data.frame(apply(x, 2, min))
colnames(raw.gwas) <- "pvalue"

raw.gwas <- as.data.frame(raw.gwas[order(raw.gwas$pvalue), ])

raw.gwas <- read.csv("raw.gwas.csv")


#OMNI GWAS result
setwd("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Data for sorghum/sorghum/Results")
omni.gwas <- read.csv("pvalue.combine.OMNIBUS.csv")
x <- as.data.frame(apply(x, 1, min))

