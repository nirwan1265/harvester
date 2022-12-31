## Loading Numerical hapmap genotype file
#Numerical file should contain the only the taxa present in the phenotype.
#this is done by filtering with taxa in tassel
#Put the taxa name on the bar in tassel and filter

#Save as numerical in TASSEL
#IMPORTANT: This can be done in TASSEL. Remove <marker> and <numerical> text from the transposed txt file before loading 
#pre processing step
#Remove the first Element <Marker>
# $ cut -f2- chr1.txt > chr1.num.txt
#Remove the first row 
# $ sed -e '1d' < chr1.num.txt > chr1.1num.txt
#Separate the file in header and body
#awk 'FNR <=1' sb_snpsDryad_sept2013_filter.c01_numerical.txt > header1.txt
#sed '1d' sb_snpsDryad_sept2013_filter.c01_numerical.txt > body1.txt
#Script to change NA to 9, 0 to 2, 1 to 0 and 0.5 to 1 and then 2.5 to 1 because it replaces 0.5 as 2.5 
# $ perl -pi -e 's/NA/9/g' body1.txt
#combine the header and body
# cat header1.txt body1.txt > originalfile.txt


##Reading the genotype files


##Reading the genotype files
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Lasky.hapmap/raw/africa.filtered")
#edit the files like the one from geno
#ch <- c("ch2L","ch2R","ch3L","ch3R","chX")
ch = c(1,2,3,4,5,6,7,8,9,10)
for(i in ch){
  assign(paste0("geno_",i), read.table(file = paste0(i,".pheno.MAF10.miss20.filtered.txt"), header = TRUE, sep = "\t"))
}


 ##Saving the genotype file as  RDS
# j <- 1
# for(i in paste0("geno", sprintf("%02d", 1:10))){
#   d = get(i)
#   saveRDS(d, paste0("geno_f",sprintf("%02d" , j),".RDS"))
#   assign(i,d)
#   j <- j+1
# }

## Loading the genotype file 

