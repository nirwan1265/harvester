install.packages("smartsnp")
library("smartsnp")

# Path to example genotype matrix "dataSNP"
pathToGenoFile = system.file("extdata", "dataSNP", package = "smartsnp")
pathToGenoFile
x <- read.table(file = pathToGenoFile)

# Example 1: modern samples
#assign 50 samples to each of two groups and colors
my_groups <- c(rep("A", 50), rep("B", 50)); cols = c("red", "blue")
my_groups

sample.id <- as.matrix(sample.id)
sample.id

x <- "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Data for sorghum/sorghum/Lasky.hapmap/eigenstratG.eg.txt"
#run PCA with truncated SVD (PCA 1 x PCA 2)
pcaR1 <- smart_pca(snp_data = x, sample_group = sample.id)


x <- pcaR1$pca.eigenvalues # extract eigenvalues
x

pcaR1$pca.snp_loadings # extract principal coefficients (SNP loadings)
pcaR1$pca.sample_coordinates # extract principal components (sample position in PCA space)

#plot PCA
plot(pcaR1$pca.sample_coordinates[,c("PC1","PC2")], cex = 2,
     pch = 19, col = cols[as.factor(my_groups)], main = "genotype smartpca")
legend("topleft", legend = levels(as.factor(my_groups)), cex =1,
       pch = 19, col = cols, text.col = cols)
# Example 2: modern and ancient samples (ancient samples projected onto modern PCA space)
#assign samples 1st to 10th per group to ancient
my_ancient <- c(1:10, 51:60)
#run PCA with truncated SVD (PCA 1 x PCA 2)
pcaR2 <- smart_pca(snp_data = pathToGenoFile, sample_group = my_groups, sample_project = my_ancient)
pcaR2$pca.eigenvalues # extract eigenvalues
pcaR2$pca.snp_loadings # extract principal coefficients (SNP loading)
pcaR2$pca.sample_coordinates # extract principal components (sample position in PCA space)
#assign samples to groups (A, ancient, B) and colors
my_groups[my_ancient] <- "ancient"; cols = c("red", "black", "blue")
#plot PCA
plot(pcaR2$pca.sample_coordinates[,c("PC1","PC2")],
     cex = 2, col = cols[as.factor(my_groups)], pch = 19, main = "genotype smartpca")
legend("topleft", legend = levels(as.factor(my_groups)), cex = 1,
       pch = 19, col = cols, text.col = cols)