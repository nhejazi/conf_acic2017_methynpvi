# get sites with raw p-values falling below cutoff
test_outcome <- 5  # choose this neurodevelopmental outcome arbitrarily
cutoff <- 0.01
tt_psiq7y <- tt_neuro_full[[test_outcome]]
targetSites <- which(tt_psiq7y$P.Value < cutoff)

# scaling of outcome based on Targeted Learning, section 7.2.2, pp. 124-125
outcome <- designMats_full[[test_outcome]][, 2]

# set up input methylation 450k data
measures450k <- as.data.frame(mcols(cpg_gr_mval))

# set up design matrix with missing cases for function input
designMatrix <- subset(batch,
                       select = c(outcomes[test_outcome], confounders, "batch"))

# reduce CpG sites based on clusters...
cpg_gr_mval <- cpg_gr
mcols(cpg_gr_mval) <- as.data.frame(cbind(clusters, mvals))
