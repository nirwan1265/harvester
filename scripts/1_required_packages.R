# Package names
packages <- c("tidyverse","ggplot2", "Rsamtools","GenomicAlignments","rtracklayer","GenomicRanges","AnnotationHub","knitr","gtools","data.table","stringi","GBJ","metap","multtest","Hmisc","devtools","SNPRelate","gdsfmt","dplyr","vcfR","tidyr","AssocTests","SKAT","NCmisc","ACAT","PANTHER.db","UniProt.ws","ape","raster","sp","rgdal","rworldmap","janitor")

# Install packages not yet installed
# installed_packages <- packages %in% rownames(installed.packages())
# if (any(installed_packages == FALSE)) {
#  install.packages(packages[!installed_packages])
# }
# if (any(installed_packages == FALSE)) {
#   BiocManager::install(packages[!installed_packages])
# }
# devtools::install_github("yaowuliu/ACAT")
#install.packages("/Users/nirwantandukar/Documents/Sorghum root rnaseq data_low phosphorus/org.Sbicolor.eg.db", repos=NULL, type="source")


# Packages loading
invisible(lapply(packages, library, character.only = TRUE))


