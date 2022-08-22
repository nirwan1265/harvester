#Loading query SNPs for GWAS from RDS file

##Location in the server: /rsstu/users/r/rrellan/sara/SorghumGEA/results/GLM_20220224 - GLM
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/GWAS.results/MLM")
##Reading the gwas RDS files
for(i in sprintf("%02d", 1:10)){
  assign(paste0("query.snp.gwas", i) , readRDS(file = paste0("mlm_sol_VL_",i,"_20220307_19_13.RDS")))
}

##Subsetting the required columns for analysis
for(i in sprintf("%02d", 1:10)){
  assign(paste0("query.snp.gwas",i), data.frame(get(paste0("query.snp.gwas",i))[3]))
  assign(paste0("query.snp.gwas",i), get(paste0("query.snp.gwas",i))[,c(2,3,4,6,7)])
}


