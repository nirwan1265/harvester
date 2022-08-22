#Required Packages
library(clusterProfiler)
library(enrichplot)
library(ggnewscale)
library(ggupset)
library(AnnotationDbi)
library(AnnotationHub)
library(GenomeInfoDb)
library(org.Sbicolor.eg.db)

# Loading gene-based pvalue data
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/pvalues.combination")
filtered_genes_OMNI <- read.csv("filtered_genes_OMNI.csv")
filtered_genes_CCT <- read.csv("filtered_genes_CCT.csv")
filtered_genes_GBJ <- read.csv("filtered_genes_GBJ.csv")
filtered_genes_SKAT <- read.csv("filtered_genes_SKAT.csv")
filtered_genes_minP <- read.csv("filtered_genes_minP.csv")
filtered_genes_GHC <- read.csv("filtered_genes_GHC.csv")


# Loading raw GWAS result
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Lasky.hapmap/R_saved")

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


#Combining total pvalues in a dataframe
raw.genes <- NULL
for (i in paste0("gwas",sprintf("%02d", 1:10),".pvalue")){
  d = get(i)
  raw.genes <- rbind(raw.genes,t(d[1,]))
  assign(i,d)
}
gene.names <- rownames(raw.genes)
raw.genes <-as.data.frame(as.numeric(raw.genes))
rownames(raw.genes) <- gene.names
colnames(raw.genes) <- "gene"
write.csv(gene.names,"gene.names.csv")


#Converting to Uniprot
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Sorghum.annotation/Phytozome/PhytozomeV12/Sbicolor/annotation")
synonym <- read.delim("Sbicolor_454_v3.1.1.synonym.txt", header = F)
colnames(synonym)[1] <- "gene.names"
#synonym <- synonym[!duplicated(synonym$V2), ]

#Converting locus to transcript
trans <- read.delim("Sbicolor_454_v3.1.1.locus_transcript_name_map.txt", header = T)
colnames(trans)[1] <- "gene.names"
#trans <- trans[!duplicated(trans$gene.names), ]
gene.names <- as.data.frame(gene.names)
colnames(gene.names) <- "gene.names"

trans2 <- inner_join(trans, gene.names)
colnames(trans2)[3] <- "gene.names"
colnames(trans2)[1] <- "gene.names."
transformed.gene <- inner_join(trans2, synonym)
transformed.gene <- transformed.gene[!duplicated(transformed.gene$V2), ]
transformed.gene <- transformed.gene$V2
write.csv(transformed.gene,"transformed.gene.csv")

AnnotationDbi::select(org.Sbicolor.eg.db, keys=transformed.gene, columns='ENTREZID', keytype='GENENAME')

#GO analysis

#Sorting the data
## feature 1: numeric vector
geneList = rnaseq_analysis[,2]
geneList = raw.genes[,1]

## feature 2: named vector
names(geneList) = as.character(rnaseq_analysis[,1])
geneList
names(geneList) = as.character(rownames(raw.genes))
geneList

## feature 3: decreasing order
geneList = sort(geneList, decreasing = TRUE)
head(geneList)

gene <-names(geneList)
head(gene)
 

#GO over-representation analysis
ego_BP <- enrichGO(gene = gene,
                   OrgDb         = org.Sbicolor.eg.db,
                   keyType       = 'ENTREZID',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
?enrichGO()
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

#Dotplot
dotplot(ego_BP, showCategory=100)
dotplot(ego_MF, showCategory=100)
dotplot(ego_CC, showCategory=100)


#Gene Concept Network:
p1_BP <- cnetplot(ego_BP, foldChange = geneList)
p3_BP <- cnetplot(ego_BP, foldChange = geneList, circular=TRUE, colorEdge=TRUE)
cowplot::plot_grid(p1_BP,p3_BP, ncol=2, labels=LETTERS[1:2], rel_widths = c(.8,.8,1.2))


p1_MF <- cnetplot(ego_MF, foldChange = geneList)
p3_MF <- cnetplot(ego_MF, foldChange = geneList, circular=TRUE, colorEdge=TRUE)
cowplot::plot_grid(p1_MF,p3_MF, ncol=2, labels=LETTERS[1:2], rel_widths = c(.8,.8,1.2))


p1_CC <- cnetplot(ego_CC, foldChange = geneList)
p3_CC <- cnetplot(ego_CC, foldChange = geneList, circular=TRUE, colorEdge=TRUE)
cowplot::plot_grid(p1_CC,p3_CC, ncol=2, labels=LETTERS[1:2], rel_widths = c(.8,.8,1.2))

#GO plot
goplot(ego_BP)
goplot(ego_MF)
goplot(ego_CC)




######################################################################################################################################################
#KEGG
######################################################################################################################################################
search_kegg_organism('ghum', by='scientific_name')
#kegg_code= sbi

sbi <- search_kegg_organism('Sorghum bicolo', by='scientific_name')
dim(sbi)
head(sbi)

kk <- enrichKEGG(gene = gene,
                 organism = 'sbi',
                 pvalueCutoff = 0.05)
head(kk)
dotplot(kk)

#GSEA
kk2 <- gseKEGG(geneList = geneList,
               organism = 'sbi',
               minGSSize = 120,
               pvalueCutoff = 0.05,
               verbose = F)
dotplot(kk2)


#KEGG module over-representation analysis
mkk <- enrichMKEGG(gene=gene,
                   organism = 'sbi',
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)
dotplot(mkk)

mkk2 <- gseMKEGG(geneList = geneList,
                 organism = 'sbi',
                 pvalueCutoff = 1)
dotplot(mkk2)


#Visualize enriched KEGG
browseKEGG(kk,'sbi00940')

library(pathview)
sbi00940 <- pathview(gene.data = geneList,
                     pathway.id = "sbi00940",
                     species = "sbi",
                     limit = list(gene=max(abs(geneList)), cpd=1))


#Reactome pathway over-representation analysis
library(remotes)
remotes::install_github("guokai8/Enrichr")
library(EnrichR)
gse <- gsea(rnaseq_analysis, SbicolorSC187_694_v1.1.annotation_info.txt, annot.info = NULL, minSize = 15, maxSize = 500, padj.method = "BH" )
gsea()

#GSEA
edo <- gseGO(geneList, 
             OrgDb  = org.Sbicolor.eg.db,
             keyType       = 'ENTREZID',
             ont           = "BP",
             pAdjustMethod = "BH",
             pvalueCutoff  = 0.01
)
edox <- setReadable(edo, 'org.Sbicolor.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=geneList)
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList)
p3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE) 
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))


p1 <- cnetplot(edox, node_label="category", 
               cex_label_category = 1.2) 
p2 <- cnetplot(edox, node_label="gene", 
               cex_label_gene = 0.8) 
p3 <- cnetplot(edox, node_label="all") 
p4 <- cnetplot(edox, node_label="none", 
               color_category='firebrick', 
               color_gene='steelblue') 
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])


p1 <- heatplot(edox, showCategory=5)
p2 <- heatplot(edox, foldChange=geneList, showCategory=5)
cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])


edox2 <- pairwise_termsim(edox)
p1 <- treeplot(edox2)
p2 <- treeplot(edox2, hclust_method = "average")
aplot::plot_list(p1, p2, tag_levels='A')

ridgeplot(edo)

gseaplot2(edo, geneSetID = 1, title = edo$Description[1])
gseaplot2(edo, geneSetID = 1:3)


library(ggplot2)
library(cowplot)
pp <- lapply(1:3, function(i) {
  anno <- edo[i, c("NES", "pvalue", "p.adjust")]
  lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")
  
  gsearank(edo, i, edo[i, 2]) + xlab(NULL) +ylab(NULL) +
    annotate("text", 10000, edo[i, "enrichmentScore"] * .75, label = lab, hjust=0, vjust=0)
})
plot_grid(plotlist=pp, ncol=1)





edo <- gseGO(geneList, 
             OrgDb  = org.Sbicolor.eg.db,
             keyType       = 'ENTREZID',
             ont           = "MF",
             pAdjustMethod = "BH",
             pvalueCutoff  = 0.01
)
edox <- setReadable(edo, 'org.Sbicolor.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=geneList)
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList)
p3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE) 
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))



write.csv(raw.gwas.gene,"raw.gwas.gene.csv")
write.csv(omnibus.gwas.gene,"omnibus.gwas.gene.csv")
system("pwd")

write.csv(gwas01.fstat,"gwas01.fstat.csv")
write.csv(gwas01.Marker,"gwas01.marker.csv")
write.csv(gwas01.pvalue,"gwas01.pvalue.csv")
write.csv(tab.pc01,"tab.pc01.csv")


