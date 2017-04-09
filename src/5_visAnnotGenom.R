# adapted from a tutorial on visualizing DNA methylation data
# https://htmlpreview.github.io/?https://github.com/PeteHaitch/tutorial.450k/blob/master/inst/doc/methylation450k.html#plotting-dmps

library(Gviz)
genome <- "hg19"
FDRcutoff = 0.05

# for visualizations, we'll use output from methytmle (but not run the function)
load(here("data", "methytmleResults2017Feb09.RData"))

# create subset of CpG sites to visualize
cpgsFDRtmle <- methytmle[which(methytmle$pvalFDR < FDRcutoff), ]
tmleCpGs <- cpg_gr_mval[which(methytmle$pvalFDR < FDRcutoff)]
mcols(tmleCpGs) <- cpgsFDRtmle

# construct a plot for each Chromosome...
targetChr <- 1
tmleCpGsToPlotGR <- tmleCpGs[seqnames(tmleCpGs) == targetChr]

chrom <- paste0("chr", as.character(seqnames(tmleCpGsToPlotGR)))
start <- start(ranges(tmleCpGsToPlotGR))
end <- end(ranges(tmleCpGsToPlotGR))
minbase <- min(start) - 0.25
maxbase <- max(end) + 0.25
pal <- c("#E41A1C", "#377EB8")

# Start building the tracks
iTrack <- IdeogramTrack(genome = genome,
                        chromosome = chrom,
                        name = "")
gTrack <- GenomeAxisTrack(col = "black",
                          cex = 1,
                          name = "",
                          fontcolor = "black")
# NOTE: This track takes a little while to create
rTrack <- UcscTrack(genome = genome,
                    chromosome = chrom,
                    track = "refGene",
                    from = min(start),
                    to = max(end),
                    trackType = "GeneRegionTrack",
                    rstarts = "exonStarts",
                    rends = "exonEnds",
                    gene = "name",
                    symbol = "name2",
                    transcript = "name",
                    strand = "strand",
                    fill = "darkblue",
                    stacking = "squish",
                    name = "RefSeq",
                    showId = TRUE,
                    geneSymbol = TRUE)
# methylation data track
gr <- tmleCpGsToPlotGR
mcols(gr) <- NULL
gr$ATE <- mcols(tmleCpGsToPlotGR)$Param
methTrack <- DataTrack(range = gr3,
                       genome = genome,
                       chromosome = chrom,
                       col = pal,
                       data = as.matrix(t(as.matrix(mcols(gr3)))),
                       ylim = c(min(gr3$lowerCI) - 1, max(gr3$upperCI) + 1),
                       type = "p",
                       name = "Local Average Tx Effect \n(DNA Methylation)",
                       background.panel = "white",
                       legend = TRUE,
                       cex.title = 0.8,
                       cex.axis = 0.8,
                       cex.legend = 0.8)

# Finally, plot the tracks
tracks <- list(iTrack, gTrack, methTrack, rTrack)
sizes <- c(2, 2, 5, 3) # set up the relative sizes of the tracks
pdf(paste0(normalizePath(here("..", "graphs")), "/plotChr", targetChr, ".pdf"),
    width = 10, height = 10)
plotTracks(tracks,
           from = minbase,
           to = maxbase,
           showTitle = TRUE,
           add53 = TRUE,
           add35 = TRUE,
           grid = TRUE,
           lty.grid = 3,
           sizes = sizes,
           length(tracks))
dev.off()
