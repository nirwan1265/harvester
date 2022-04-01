# `harvester`
Combine candidate gene evidence from multiple sources

Reference genomes

```{sh}
# in csh
tree /rsstu/users/r/rrellan/sara/ref

# /rsstu/users/r/rrellan/sara/ref
# ├── B73_RefGen_v2.fa
# ├── B73_RefGen_v2.fa.fai
# ├── Zm-B73-REFERENCE-GRAMENE-4.0.fa
# ├── Zm-B73-REFERENCE-GRAMENE-4.0.fa.fai
# ├── Zm-B73-REFERENCE-NAM-5.0.fa
# └── Zm-B73-REFERENCE-NAM-5.0.fa.fai
```

## Data to liftover into v5

As there is no v2 to v5 chainfile, that liftover must be done in two steps v2 -> v4 and v4 -> v5

```{sh}
# in csh
tree /rsstu/users/r/rrellan/sara/ZeaGEA/results

# /rsstu/users/r/rrellan/sara/ZeaGEA/results
# ├── glm_20181019
# │   ├── GLM_Q28_site                 # GWAS pvalues AGPv2 corrected with K=28 subpopulations in ADMIXTURE
# │   └── GLM_site                     # GWAS pvalues AGPv2 no correction for kinship or population structure
# └── pcadapt
#     ├── AGPv2_to_B73_RefGen_v4.chain
#     ├── B73_RefGen_v2_Chr.bed
#     ├── pcadapt_corrected.Rimage     # PCAdapt pvalues AGPv2
#     └── PCAdapt_v4.RDS               # PCAdapt pvalues AGPv4
```

`B73_RefGen_v2_Chr.bed` BED file is built like:

```{sh}
# in csh
# activate conda environment for samtools
conda activate /usr/local/usrapps/maize/sorghum/conda/envs/snakemake-tutorial

set ref=/rsstu/users/r/rrellan/sara/ref
cd $ref 

wget https://ftp.maizegdb.org/MaizeGDB/FTP/B73_RefGen_v2/B73_RefGen_v2.fa.gz 
gunzip B73_RefGen_v2.fa.gz
samtools faidx B73_RefGen_v2.fa

set pcadapt=/rsstu/users/r/rrellan/sara/ZeaGEA/results/pcadapt

cd $pcadapt

head -n 10 $ref/B73_RefGen_v2.fa.fai | awk 'BEGIN {FS="\t"}; {gsub("Chr", "", $1); gsub(":.*", "", $1); print $1 FS "0" FS $2 FS "ws"}' | sort -k1,1n > B73_RefGen_v2_Chr.bed

```

Chain files [(see specs)](https://genome.ucsc.edu/goldenPath/help/chain.html)  downloaded from [the MaizeGDB ftp site](http://ftp.gramene.org/CURRENT_RELEASE/assembly_chain/zea_mays/) won't be read by R unless spaces are replaced by tabs:

```{sh}
# in csh
wget "http://ftp.gramene.org/CURRENT_RELEASE/assembly_chain/zea_mays/AGPv3_to_B73_RefGen_v4.chain.gz"
gunzip AGPv3_to_B73_RefGen_v4.chain.gz
perl -i -pe ' if ( $_ !~ /chain/) {s/ +/\t/g}' AGPv3_to_B73_RefGen_v4.chain
```

## Script to change NA to 0 for numerical genotype file
```{sh}
# in zsh
perl -pi -e 's/NA/0/g' numerical.genotype.chr<chr_no.>.txt
```

## Other possible data for combining:

  - PBE scores (these are not p-values)
  
  - Differential expression in +P -P experiments



