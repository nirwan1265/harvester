library(GBJ)
# Load the genotype data at FGFR2 SNPs for 91 Great Britain (GBR) subjects from the 1000 Genomes Project (publically available).
data(FGFR2)

# Load PCs for these same 91 subjects (calculated, for example, with EIGENSTRAT)
data(gbr_pcs)
gbr_pcs
typeof(gbr_pcs)
class(gbr_pcs)
# Suppose we were given the following 64 test statistics for the FGFR2 SNPs (must be the same SNPs that
# are in our genotype matrix!)
FGFR2_stats <- rnorm(n=64)
FGFR2
typeof(FGFR2)
class(FGFR2)


# Estimate correlation matrix for summary statistics
FGFR2_cor_mat <- estimate_ss_cor(ref_pcs=gbr_pcs, ref_genotypes=FGFR2, link_function='logit')
gbr_pcs
FGFR2

typeof(FGFR2)
class(FGFR2)


# Run GBJ
GBJ(test_stats=FGFR2_stats, cor_mat=FGFR2_cor_mat)