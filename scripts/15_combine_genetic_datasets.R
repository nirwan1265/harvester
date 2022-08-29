# Loading Data
## ONMI data
omni_uniprot <- read.csv("OMNI_all_filtered.for.uniprot.csv")
head(omni_uniprot)
omni <- omni_uniprot[,c(1,3)]
omni <- unique(omni)
head(omni)

## RNAseq data
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/RNAseq_lowP")
rnaseq_unirprot <- read.csv("lowP_sorghum_rnaseq_diff.csv")
head(rnaseq_unirprot)
rnaseq <- rnaseq_unirprot[,c(3,8)]
head(rnaseq)
colnames(rnaseq) <- c("pvalue", "Gene")
rnaseq[rnaseq == 0] <- min(rnaseq$pvalue[rnaseq$pvalue > 0])/2
rnaseq <- rnaseq[!duplicated(rnaseq$Gene), ]




#Combining datasets using Cauchy Combination Test:
## Data prep
omni_rnaseq <- inner_join(omni,rnaseq, by = "Gene")
head(omni_rnaseq)
colnames(omni_rnaseq) <- c("GeneName", "pvalue.omni","pvalue.rnaseq")


## Using function
omni_rnaseq_cct <- cbind(omni_rnaseq, (apply(omni_rnaseq[,c(2,3)],1,acat)))
colnames(omni_rnaseq_cct)[4] <- "pvalue.CCT"
omni_rnaseq_cct_filter <- omni_rnaseq_cct[which(omni_rnaseq_cct$pvalue.CCT < 0.05), ]
head(omni_rnaseq_cct_filter)
write.csv(omni_rnaseq_cct_filter,"omni_rnaseq_cct_filter.csv")


# Loading converted gene files for CCT_gwas_rnaseq_gene.conversion_filtered.csv
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/pvalues.combination")
gwas_rnaseq <- read.csv("CCT_gwas_rnaseq_gene.conversion_filtered.csv")



#GO analysis

#Sorting the data
## feature 1: numeric vector
#geneList = rnaseq_analysis[,2]
geneList = gwas_rnaseq[,4]

## feature 2: named vector
#names(geneList) = as.character(rnaseq_analysis[,1])
#geneList
names(geneList) = as.character(gwas_rnaseq[,3])
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
go.bp.res %>% 
  ggplot(aes(Description,Count))+
  geom_col() +
  coord_flip()+
  labs(title="Enriched Biological Process using CCT for eGWAS for Phosphorus and root RNAseq on low Phosphours",
       x="Gene Count", y= "Description")+
  geom_text(aes(label = round(Count, 1)), nudge_y= -1, color="white")

go.mf.res <- as.data.frame(cbind(ego_MF@result[["Description"]],ego_MF@result[["Count"]],ego_MF@result[["pvalue"]]))
colnames(go.mf.res) <- c("Description","Count","pvalue")
go.mf.res[,3] <- sapply(go.mf.res[,3], as.numeric)
go.mf.res <- go.mf.res[which(go.mf.res$pvalue < 0.05), ]
go.mf.res[,2] <- sapply(go.mf.res[,2], as.double)
go.mf.res %>% 
  ggplot(aes(Description,Count))+
  geom_col() +
  coord_flip()+
  labs(title="Enriched Molecular Functions using CCT for eGWAS for Phosphorus and root RNAseq on low Phosphours",
       x="Gene Count", y= "Description")+
  geom_text(aes(label = round(Count, 1)), nudge_y= -1, color="white")

go.cc.res <- as.data.frame(cbind(ego_CC@result[["Description"]],ego_CC@result[["Count"]],ego_CC@result[["pvalue"]]))
colnames(go.cc.res) <- c("Description","Count","pvalue")
go.cc.res[,3] <- sapply(go.cc.res[,3], as.numeric)
go.cc.res <- go.cc.res[which(go.cc.res$pvalue < 0.05), ]
go.cc.res[,2] <- sapply(go.cc.res[,2], as.double)
go.cc.res %>% 
  ggplot(aes(Description,Count))+
  geom_col() +
  coord_flip()+
  labs(title="Enriched Cellular Processes using CCT for eGWAS for Phosphorus and root RNAseq on low Phosphours",
       x="Gene Count", y= "Description")+
  geom_text(aes(label = round(Count, 1)), nudge_y= -1, color="white")
