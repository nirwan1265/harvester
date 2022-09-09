# Loading files:

#Z stat
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/R_saved/filtered.gwas")
for(i in sprintf("%02d", 1:10)){
  assign(paste0("gwas",i,".zstat"), readRDS(paste0("gwas",i,".zstat.RDS")))
}

# pvalue
for(i in sprintf("%02d", 1:10)){
  assign(paste0("gwas",i,".pvalue"), readRDS(paste0("gwas",i,".pvalue.RDS")))
}

# Marker
for(i in sprintf("%02d", 1:10)){
  assign(paste0("gwas",i,".Marker"), readRDS(paste0("gwas",i,".Marker.RDS")))
}

# Genotype
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/R_saved/genotype_filtered.phenotype")
for(i in sprintf("%02d", 1:10)){
  assign(paste0("geno_f",i), readRDS(paste0("geno_f",i,".RDS")))
}


# PCA
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/R_saved/pca")
for(i in sprintf("%02d", 1:10)){
  assign(paste0("tab",i), readRDS(paste0("tab",i,".RDS")))
}


# Phenotype
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/R_saved/phenotype")
pheno <- readRDS("pheno.RDS")
