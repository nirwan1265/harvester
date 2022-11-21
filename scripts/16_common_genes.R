# Working directory
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/pvalues.combination/MAGMA/sorghum/combined.pvalues.snpwise")



#Loading files
folder <- getwd()
file_list <- list.files(path = folder)

for(i in 1:length(file_list)){
  assign(file_list[i], read.delim(file_list[i], sep = ""))
}
rm(`0`)



tot <- read.delim("tot_sorghum_LMM.magma.multi.snpwise.genes.out", sep = "")

#Loading all the files together in one common variable name

files <- vroom(fs::dir_ls(glob = "*genes.out"))
files
data <- vroom(files)

head(data)
