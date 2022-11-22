# Working directory
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/pvalues.combination/MAGMA/sorghum/combined.pvalues.snpwise")

#Loading files
folder <- getwd()
file_list <- list.files(path = folder, pattern = ".txt")
file_list
for(i in 1:length(file_list)){
  assign(file_list[i], read.delim(file_list[i], sep = ""))
}
for(i in 1:length(file_list)){
  d = get(file_list[i])
  print(d[1:4,1:4])
  d = d[which(d[9] < 0.05), ]
  d = d[,c(1,9)]
  #assign(file_list[i])
  assign(file_list[i],d)
}


# Making a list of all the phenotypes:
list_pheno <- list(apa = apa_sorghum_LMM.magma.multi.snpwise.genes.out.txt$GENE,
                   ext = ext_P20_sorghum_LMM.magma.multi.snpwise.genes.out.txt$GENE, 
                   lab = PBR1_sorghum_LMM.magma.multi.snpwise.genes.out.txt$GENE, 
                   NPlim = NPlim_sorghum_LMM.magma.multi.snpwise.genes.out.txt$GENE, 
                   occ = occ_sorghum_LMM.magma.multi.snpwise.genes.out.txt$GENE, 
                   org = org_sorghum_LMM.magma.multi.snpwise.genes.out.txt$GENE, 
                   PBR1 = PBR1_sorghum_LMM.magma.multi.snpwise.genes.out.txt$GENE, 
                   PBR2 = PBR2_sorghum_LMM.magma.multi.snpwise.genes.out.txt$GENE, 
                   PHO1 = PHO1_sorghum_LMM.magma.multi.snpwise.genes.out.txt$GENE, 
                   PHO2 = PHO2_sorghum_LMM.magma.multi.snpwise.genes.out.txt$GENE,
                   PMEH1 = PMEH1_sorghum_LMM.magma.multi.snpwise.genes.out.txt$GENE, 
                   PMEH2 = PMEH2_sorghum_LMM.magma.multi.snpwise.genes.out.txt$GENE, 
                   PNZ1 = PNZ1_sorghum_LMM.magma.multi.snpwise.genes.out.txt$GENE, 
                   PNZ2 = PNZ2_sorghum_LMM.magma.multi.snpwise.genes.out.txt$GENE, 
                   POL1 = POL1_sorghum_LMM.magma.multi.snpwise.genes.out.txt$GENE, 
                   POL2 = POL2_sorghum_LMM.magma.multi.snpwise.genes.out.txt$GENE, 
                   sec = sec_sorghum_LMM.magma.multi.snpwise.genes.out.txt$GENE, 
                   sol_Hi = sol_Hi_sorghum_LMM.magma.multi.snpwise.genes.out.txt$GENE, 
                   sol_Lo = sol_Lo_sorghum_LMM.magma.multi.snpwise.genes.out.txt$GENE, 
                   sol_Mo = sol_Mo_sorghum_LMM.magma.multi.snpwise.genes.out.txt$GENE, 
                   sol = sol_sorghum_LMM.magma.multi.snpwise.genes.out.txt$GENE, 
                   stp10 = stp10_sorghum_LMM.magma.multi.snpwise.genes.out.txt$GENE, 
                   stp10 = stp10_sorghum_LMM.magma.multi.snpwise.genes.out.txt$GENE, 
                   stp20 = stp20_sorghum_LMM.magma.multi.snpwise.genes.out.txt$GENE, 
                   stp30 = stp30_sorghum_LMM.magma.multi.snpwise.genes.out.txt$GENE, 
                   stp100 = stp100_sorghum_LMM.magma.multi.snpwise.genes.out.txt$GENE, 
                   tot = tot_sorghum_LMM.magma.multi.snpwise.genes.out.txt$GENE, 
                   TP1 = TP1_sorghum_LMM.magma.multi.snpwise.genes.out.txt$GENE, 
                   TP2 = TP2_sorghum_LMM.magma.multi.snpwise.genes.out.txt$GENE) 

x = list("apa", "ext", "TP1")

x = list(expand.grid(c("apa","ext","lab","NPlim","occ","TP1","TP2"), c("apa","ext","lab","TP1","TP2"),"lab", "TP1","TP2"))
outer(c("apa","ext","lab","TP1","TP2"), c("apa","ext","lab","TP1","TP2"), "paste")

x = c(c("apa","ext","lab","org","PBR1","PBR2","TP1","TP2"), c("apa","ext","lab","org","PBR1","PBR2","TP1","TP2"))
x = unique(x)
combn(x,7)

x
# Finding Common genes in sets of two: 
nms <- combn( names(list_pheno) , 2 , FUN = paste0 , collapse = "" , simplify = FALSE )
nms

ll <- combn( list_pheno , 2 , simplify = FALSE )
ll
out <- lapply( ll , function(x)  intersect( x[[1]] , x[[2]] ) ) 
out_length <- lapply( ll , function(x)  length(intersect( x[[1]] , x[[2]] ) ))
length_commonphenot <- setNames(out_length, nms)
length_commonphenot <- lapply(length_commonphenot, sort)
