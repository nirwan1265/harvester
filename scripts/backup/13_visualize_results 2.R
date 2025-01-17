#Filtering genes
#Raw GWAS result
setwd("/Users/nirwan/Library/Mobile Documents/com~apple~CloudDocs/Data for sorghum/sorghum/Results")
x <- gwas01.pvalue
x[] <- lapply(x, function(x) as.numeric(replace(x, is.na(x), 1)))

raw.gwas <- as.data.frame(apply(x, 2, min))
colnames(raw.gwas)[1] <- "pvalue"
raw.gwas$gene <- rownames(raw.gwas)
raw.gwas <- as.data.frame(raw.gwas[order(raw.gwas$pvalue), ])
raw.gwas <- raw.gwas[raw.gwas$pvalue < 0.05, ]
raw.gwas.gene <- raw.gwas$gene




#OMNI GWAS result
setwd("/Users/nirwan/Library/Mobile Documents/com~apple~CloudDocs/Data for sorghum/sorghum/Results")
omni.gwas <- read.csv("pvalue.combine.OMNIBUS.csv")
x <- omni.gwas
x <- as.data.frame(apply(x, 1, FUN = min))
x[] <- lapply(x, function(x) as.numeric(replace(x, is.na(x), 1)))
omnibus.gwas <- cbind(omni.gwas,x)
omnibus.gwas <- omnibus.gwas[,c(1,6)]
colnames(omnibus.gwas) <- c("gene","pvalue")
omnibus.gwas <- as.data.frame(omnibus.gwas[order(omnibus.gwas$pvalue), ])

omnibus.gwas <- omnibus.gwas[omnibus.gwas$pvalue < 0.05, ]
omnibus.gwas.gene <-omnibus.gwas$gene

library(dplyr)
x <- left_join(omnibus.gwas.gene, raw.gwas.gene)

#Venn diagram
if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")
library(ggvenn)

x <- list(
  pvalue.combination = omnibus.gwas.gene,
  RAW.GWAS = raw.gwas.gene
)

ggvenn(
  x,
  fill_color = c("#CD534CFF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
)

#GO analysis
library(clusterProfiler)
library(enrichplot)
library(ggnewscale)
library(ggupset)
library(AnnotationDbi)
library(AnnotationHub)
library(GenomeInfoDb)
library(org.Sbicolor.eg.db)



#Sorting the data
## feature 1: numeric vector
geneList = rnaseq_analysis[,2]
geneList

## feature 2: named vector
names(geneList) = as.character(rnaseq_analysis[,1])
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
