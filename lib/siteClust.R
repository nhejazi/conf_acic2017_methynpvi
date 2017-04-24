library(bumphunter)

clusterSites <- function(granges, ) {
  if(dim(mcols(granges))[2] != 0) {
    mcols(granges) <- NULL
  }
  clusters <- bumphunter::boundedClusterMaker(chr = seqnames(granges),
                                              pos = start(ranges(granges)),
                                              assumeSorted = FALSE,
                                              maxClusterWidth = 3000,
                                              maxGap = 300)
  return(clusters)
}
