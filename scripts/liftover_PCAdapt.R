library(liftOver)
library(dplyr)
library(GenomicRanges)
library(regioneR)

# Loading "colm" dataframe with PCAdapt pvalues
load("pcadapt_corrected.Rimage", verbose = TRUE)
PCAdapt <- colm
colm <- NULL

# SEEDs data  (Romero Navarro 2017) is in AGPv2
# This version of pcadapt is also in AGPv2
# because PCAdapt data loads into a GRanges object with warnings
# about out of range SNPS if AGPv3 is used
#
# B73_RefGen_v2_Chr.bed
# file with each chromosome start and end coordinates
# this is useful for checking the genome version of the markers
#
# REMEMBER: BED file format specifications say that bed files are indexed
# at 0 (zero indexed coordinates)!!!!!!!!!!!
#
# head -n 10 /rsstu/users/r/rrellan/sara/B73_RefGen_v2.fa.fai | awk 'BEGIN {FS="\t"}; {gsub("Chr", "", $1); gsub(":.*", "", $1); print $1 FS "0" FS $2 FS "ws"}' | sort -k1,1n > B73_RefGen_v2_Chr.bed 
#
# double check SNP coordinates!
#
# I think we must uplift all coordinates to AGPv4
#
# Indexed with samtools faidx from
# https://ftp.maizegdb.org/MaizeGDB/FTP/B73_RefGen_v2/B73_RefGen_v2.fa.gz

AGPv2 <- toGRanges("B73_RefGen_v2_Chr.bed")

seqlengths(AGPv2) <- width(AGPv2)

genome(AGPv2) <- "AGPv2"

# Loading chainfile for liftover
# source:
# http://ftp.gramene.org/CURRENT_RELEASE/assembly_chain/zea_mays/
chain_file <- "AGPv2_to_B73_RefGen_v4.chain"


ch <-import.chain(chain_file)
PCAdapt_v2 <- makeGRangesFromDataFrame(PCAdapt,
                                keep.extra.columns=TRUE,
                                ignore.strand=TRUE,
                                seqinfo = seqinfo(AGPv2),  # Sanity check
                                seqnames.field="CHR",
                                start.field="BP",
                                end.field="BP") 

PCAdapt_v4 <- liftOver(PCAdapt_v2,ch) %>% unlist()
saveRDS(PCAdapt_v4, file= "PCAdapt_v4.RDS")

#Converting to V5
chain_file2 <- "B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain"
ch2 <- import.chain(chain_file2)
PCAdapt_v5 <-  liftOver(PCAdapt_v4,ch2) %>% unlist()



