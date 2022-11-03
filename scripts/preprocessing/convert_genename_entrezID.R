#https://rdrr.io/bioc/BUSpaRse/man/tr2g_gff3.html
BiocManager::install("BUSpaRse")
library(BUSpaRse)


x <- tr2g_gff3("Sorghum_bicolor.Sorghum_bicolor_NCBIv3.54.chromosome.1.gff3", write_tr2g = FALSE, get_transcriptome = FALSE, save_filtered_gff = FALSE, gene_id = "gene_id")
