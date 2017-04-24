# FDR controlling procedure from Tuglus & van der Laan (2008))

FDR_msa <- function(pvals, totalTests, ...) {
  pvalNoTest <- rep(1, totalTests - length(pvals))
  pvalsAll <- c(pvals, pvalNoTest)
  fdr_adj <- p.adjust(pvalsAll, method = "fdr")
  fdr_out <- fdr_adj[1:length(pvals)]
  return(fdr_out)
}
