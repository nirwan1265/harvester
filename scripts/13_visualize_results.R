#Required Packages
library(clusterProfiler)
library(enrichplot)
library(ggnewscale)
library(ggupset)
library(AnnotationDbi)
library(AnnotationHub)
library(GenomeInfoDb)
library(org.Sbicolor.eg.db)
library(org.Hs.eg.db)
library(dplyr)

# Loading gene-based pvalue data
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/pvalues.combination")
filtered_genes_OMNI <- read.csv("filtered_genes_OMNI.csv")
filtered_genes_CCT <- read.csv("filtered_genes_CCT.csv")
filtered_genes_GBJ <- read.csv("filtered_genes_GBJ.csv")
filtered_genes_SKAT <- read.csv("filtered_genes_SKAT.csv")
filtered_genes_minP <- read.csv("filtered_genes_minP.csv")
filtered_genes_GHC <- read.csv("filtered_genes_GHC.csv")



# Loading raw GWAS result
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/R_saved")

gwas01.pvalue <- readRDS("gwas01.pvalue.RDS")
gwas02.pvalue <- readRDS("gwas02.pvalue.RDS")
gwas03.pvalue <- readRDS("gwas03.pvalue.RDS")
gwas04.pvalue <- readRDS("gwas04.pvalue.RDS")
gwas05.pvalue <- readRDS("gwas05.pvalue.RDS")
gwas06.pvalue <- readRDS("gwas06.pvalue.RDS")
gwas07.pvalue <- readRDS("gwas07.pvalue.RDS")
gwas08.pvalue <- readRDS("gwas08.pvalue.RDS")
gwas09.pvalue <- readRDS("gwas09.pvalue.RDS")
gwas10.pvalue <- readRDS("gwas10.pvalue.RDS")


#Combining total  Raw pvalues in a dataframe
raw.genes <- NULL
for (i in paste0("gwas",sprintf("%02d", 1:10),".pvalue")){
  d = get(i)
  raw.genes <- rbind(raw.genes,t(d[1,]))
  assign(i,d)
}
gene.names <- rownames(raw.genes)
raw.genes <-as.data.frame(as.numeric(raw.genes))
rownames(raw.genes) <- gene.names
colnames(raw.genes) <- "pvalue"
raw.genes$GeneName <- rownames(raw.genes)
#write.csv(raw.genes,"raw.gwas.csv")



# Converting Gene Names to Gene ID. Happens in two steps:
## 1.) Go to http://www.pantherdb.org/ and paste the genes, select Sorghum bicolor as the organism and click submit. 
##     Download the list of genes, with the UniProt ID and separately save them in another file. 
## 2.) Go to https://www.uniprot.org/ and change the gene name from UniProt to GeneID. Save the file. 
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Gene.conversion")
# Loading files
genename_uni <- read.csv("GeneName_UniProt.csv")
uni_geneid <- read.csv("UniProt_GeneID.csv")
genename_geneid <- inner_join(genename_uni,uni_geneid, by = "UniProt")

# Filtering out only the ones with geneid
raw.genes_geneid <- inner_join(raw.genes, genename_geneid)
rownames(raw.genes_geneid) <- raw.genes_geneid$GeneID
raw.genes_geneid <- raw.genes_geneid[,]

#GO analysis

#Sorting the data
## feature 1: numeric vector
#geneList = rnaseq_analysis[,2]
geneList = raw.genes_geneid[,1]

## feature 2: named vector
#names(geneList) = as.character(rnaseq_analysis[,1])
#geneList
names(geneList) = as.character(raw.genes_geneid[,4])
geneList

## feature 3: decreasing order
geneList = sort(geneList, decreasing = TRUE)
head(geneList)

gene <-names(geneList)[geneList < 0.05 & geneList != 0]
head(gene)




# keytype
# “ACCNUM” “ALIAS” “ENSEMBL” “ENSEMBLPROT” “ENSEMBLTRANS” “ENTREZID”
# “ENZYME” “EVIDENCE” “EVIDENCEALL” “FLYBASE” “FLYBASECG” “FLYBASEPROT”
# “GENENAME” “GO” “GOALL” “MAP” “ONTOLOGY” “ONTOLOGYALL”
# “PATH” “PMID” “REFSEQ” “SYMBOL” “UNIGENE” “UNIPROT”

#GO over-representation analysis
ego_BP <- enrichGO(gene = gene,
                   OrgDb         = org.Sbicolor.eg.db,
                   keyType       = 'ENTREZID',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)

ego_MF <- enrichGO(gene = gene,
                   OrgDb         = org.Sbicolor.eg.db,
                   keyType       = 'ENTREZID',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)

ego_CC <- enrichGO(gene = gene,
                   OrgDb         = org.Sbicolor.eg.db,
                   keyType       = 'ENTREZID',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)


# Visualize
go.bp.res <- as.data.frame(cbind(ego_BP@result[["Description"]],ego_BP@result[["Count"]],ego_BP@result[["pvalue"]]))
colnames(go.bp.res) <- c("Description","Count","pvalue")
go.bp.res[,3] <- sapply(go.bp.res[,3], as.numeric)
go.bp.res <- go.bp.res[which(go.bp.res$pvalue < 0.05), ]
go.bp.res[,2] <- sapply(go.bp.res[,2], as.double)
plot <- go.bp.res %>% 
  ggplot(aes(x = Description,y = Count,fill=-log10(pvalue)))+
  geom_col() +
  coord_flip()+
  labs(title="Raw GWAS enriched Biological Processes",
       x="Description", y= "Gene Count")+
  geom_text(aes(label = round(Count, 1)), nudge_y= -1, color="white")+
  scale_fill_viridis_b(trans='log10')+
  theme(axis.text = element_text(size = 30))+ 
  theme(axis.title = element_text(size = 25)) +
  theme(legend.text = element_text(size = 22)) +
  theme(plot.title = element_text(size = 30)) +
  theme(legend.title = element_text(size = 25)) 
ggsave("Raw GWAS enriched Biological Processes.tiff", plot, width=25, height=15, units="in", dpi=750)

  



go.mf.res <- as.data.frame(cbind(ego_MF@result[["Description"]],ego_MF@result[["Count"]],ego_MF@result[["pvalue"]]))
colnames(go.mf.res) <- c("Description","Count","pvalue")
go.mf.res[,3] <- sapply(go.mf.res[,3], as.numeric)
go.mf.res <- go.mf.res[which(go.mf.res$pvalue < 0.05), ]
go.mf.res[,2] <- sapply(go.mf.res[,2], as.double)
plot <- go.mf.res %>% 
  ggplot(aes(Description,Count,fill=-log10(pvalue)))+
  geom_col() +
  coord_flip()+
  labs(title="Raw GWAS enriched Molecular Functions",
       x="Description", y= "Gene Count")+
  geom_text(aes(label = round(Count, 1)), nudge_y= -1, color="white")+
  scale_fill_viridis_b(trans='log10')+
  theme(axis.text = element_text(size = 30))+ 
  theme(axis.title = element_text(size = 25)) +
  theme(legend.text = element_text(size = 22)) +
  theme(plot.title = element_text(size = 30)) +
  theme(legend.title = element_text(size = 25)) 
ggsave("Raw GWAS enriched Molecular Functions.tiff", plot, width=40, height=15, units="in", dpi=750)


go.cc.res <- as.data.frame(cbind(ego_CC@result[["Description"]],ego_CC@result[["Count"]],ego_CC@result[["pvalue"]]))
colnames(go.cc.res) <- c("Description","Count","pvalue")
go.cc.res[,3] <- sapply(go.cc.res[,3], as.numeric)
go.cc.res <- go.cc.res[which(go.cc.res$pvalue < 0.05), ]
go.cc.res[,2] <- sapply(go.cc.res[,2], as.double)
plot <- go.cc.res %>% 
  ggplot(aes(Description,Count,fill=-log10(pvalue)))+
  geom_col() +
  coord_flip()+
  labs(title="Raw GWAS enriched Cellular Components",
       x="Description", y= "Gene Count")+
  geom_text(aes(label = round(Count, 1)), nudge_y= -1, color="white")+
  scale_fill_viridis_b(trans='log10')+
  theme(axis.text = element_text(size = 30))+ 
  theme(axis.title = element_text(size = 25)) +
  theme(legend.text = element_text(size = 22)) +
  theme(plot.title = element_text(size = 30)) +
  theme(legend.title = element_text(size = 25)) 
ggsave("Raw GWAS enriched Cellular Components.tiff", plot, width=25, height=15, units="in", dpi=750)




################################################################################
################################################################################

# For OMNI test:
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/pvalues.combination")
## Read Data:
omni_genename_geneid <- read.csv("OMNI_filtered_genename_uniprot_geneid.csv")
colnames(omni_genename_geneid)[1] <- "Gene"
omni <- filtered_genes_OMNI[,c(1,6)]

# Changing gene name
omni_geneid <- inner_join(omni_genename_geneid, omni)

## Here a gene,C5XPA5, has several different transcripts with different geneid but same uniprot name. 
## have to put the pvalues individually for them
#Haven't done it here due to time constraint for the poster on Friday aug 23rd.


#GO analysis

#Sorting the data
## feature 1: numeric vector
#geneList = rnaseq_analysis[,2]
geneList = omni_geneid[,4]

## feature 2: named vector
#names(geneList) = as.character(rnaseq_analysis[,1])
#geneList
names(geneList) = as.character(omni_geneid[,3])
geneList

## feature 3: decreasing order
geneList = sort(geneList, decreasing = TRUE)
head(geneList)

gene <-names(geneList)[geneList < 0.05 & geneList != 0]
head(gene)




# keytype
# “ACCNUM” “ALIAS” “ENSEMBL” “ENSEMBLPROT” “ENSEMBLTRANS” “ENTREZID”
# “ENZYME” “EVIDENCE” “EVIDENCEALL” “FLYBASE” “FLYBASECG” “FLYBASEPROT”
# “GENENAME” “GO” “GOALL” “MAP” “ONTOLOGY” “ONTOLOGYALL”
# “PATH” “PMID” “REFSEQ” “SYMBOL” “UNIGENE” “UNIPROT”

#GO over-representation analysis
ego_BP <- enrichGO(gene = gene,
                   OrgDb         = org.Sbicolor.eg.db,
                   keyType       = 'ENTREZID',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)

ego_MF <- enrichGO(gene = gene,
                   OrgDb         = org.Sbicolor.eg.db,
                   keyType       = 'ENTREZID',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)

ego_CC <- enrichGO(gene = gene,
                   OrgDb         = org.Sbicolor.eg.db,
                   keyType       = 'ENTREZID',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)


# Visualize
go.bp.res <- as.data.frame(cbind(ego_BP@result[["Description"]],ego_BP@result[["Count"]],ego_BP@result[["pvalue"]]))
colnames(go.bp.res) <- c("Description","Count","pvalue")
go.bp.res[,3] <- sapply(go.bp.res[,3], as.numeric)
go.bp.res <- go.bp.res[which(go.bp.res$pvalue < 0.05), ]
go.bp.res[,2] <- sapply(go.bp.res[,2], as.double)
plot <- go.bp.res %>% 
  ggplot(aes(Description,Count,fill=-log10(pvalue)))+
  geom_col() +
  coord_flip()+
  labs(title="OMNI test enriched Biological Processes",
       x="Description", y= "Gene Count")+
  geom_text(aes(label = round(Count, 1)), nudge_y= -1, color="white")+
  scale_fill_viridis_b(trans='log10')+
  theme(axis.text = element_text(size = 30))+ 
  theme(axis.title = element_text(size = 25)) +
  theme(legend.text = element_text(size = 22)) +
  theme(plot.title = element_text(size = 30)) +
  theme(legend.title = element_text(size = 25)) 
ggsave("OMNI test enriched Biological Processes.tiff", plot, width=25, height=15, units="in", dpi=750)


go.mf.res <- as.data.frame(cbind(ego_MF@result[["Description"]],ego_MF@result[["Count"]],ego_MF@result[["pvalue"]]))
colnames(go.mf.res) <- c("Description","Count","pvalue")
go.mf.res[,3] <- sapply(go.mf.res[,3], as.numeric)
go.mf.res <- go.mf.res[which(go.mf.res$pvalue < 0.05), ]
go.mf.res[,2] <- sapply(go.mf.res[,2], as.double)
plot <- go.mf.res %>% 
  ggplot(aes(Description,Count,fill=-log10(pvalue)))+
  geom_col() +
  coord_flip()+
  labs(title="OMNI test enriched Molecular Functions",
       x="Description", y= "Gene Count")+
  geom_text(aes(label = round(Count, 1)), nudge_y= -1, color="white")+
  scale_fill_viridis_b(trans='log10')+
  theme(axis.text = element_text(size = 40))+ 
  theme(axis.title = element_text(size = 25)) +
  theme(legend.text = element_text(size = 22)) +
  theme(plot.title = element_text(size = 30)) +
  theme(legend.title = element_text(size = 25)) 
ggsave("OMNI test enriched Molecular Functions.tiff", plot, width=25, height=15, units="in", dpi=750)




go.cc.res <- as.data.frame(cbind(ego_CC@result[["Description"]],ego_CC@result[["Count"]],ego_CC@result[["pvalue"]]))
colnames(go.cc.res) <- c("Description","Count","pvalue")
go.cc.res[,3] <- sapply(go.cc.res[,3], as.numeric)
go.cc.res <- go.cc.res[which(go.cc.res$pvalue < 0.05), ]
go.cc.res[,2] <- sapply(go.cc.res[,2], as.double)
plot <- go.cc.res %>% 
  ggplot(aes(Description,Count,fill=-log10(pvalue)))+
  geom_col() +
  coord_flip()+
  labs(title="OMNI test enriched Cellular Components",
       x="Description", y= "Gene Count")+
  geom_text(aes(label = round(Count, 1)), nudge_y= -1, color="white")+
  scale_fill_viridis_b(trans='log10')+
  theme(axis.text = element_text(size = 40))+ 
  theme(axis.title = element_text(size = 25)) +
  theme(legend.text = element_text(size = 22)) +
  theme(plot.title = element_text(size = 30)) +
  theme(legend.title = element_text(size = 25)) 
ggsave("OMNI test enriched Cellular Components.tiff", plot, width=25, height=15, units="in", dpi=750)

x <- rbind(go.bp.res,go.mf.res,go.cc.res)
x %>% 
  ggplot(aes(Description,Count,fill=-log10(pvalue)))+
  geom_col() +
  coord_flip()+
  labs(title="OMNI test enriched Cellular Components",
       x="Description", y= "Gene Count")+
  geom_text(aes(label = round(Count, 1)), nudge_y= -1, color="white")+
  scale_fill_viridis_b(trans='log10')
  # theme(axis.text = element_text(size = 40))+ 
  # theme(axis.title = element_text(size = 25)) +
  # theme(legend.text = element_text(size = 22)) +
  # theme(plot.title = element_text(size = 30)) +
  # theme(legend.title = element_text(size = 25)) 
