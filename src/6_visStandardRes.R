# for visualizations, we'll use output from methytmle (but not run the function)
load(here("data", "methytmleOutput2017Feb13.RData"))
fdrCutoff = 0.05

library(ggplot2)
library(qqman)

library(wesanderson)
pal1 <- wesanderson::wes_palette("Rushmore", 100, type = "continuous")
pal2 <- wesanderson::wes_palette("Darjeeling", type = "continuous")

# created GRanges from SummarizedExperiment container output
library(SummarizedExperiment)
CpGs_gr <- rowRanges(tmleOut_se)
dmCpGs_gr <- rowRanges(tmleOut_se)[rowRanges(tmleOut_se)$pvalFDR < fdrCutoff, ]
dmCpGs_gr <- sort(dmCpGs_gr)
methytmle <- rowRanges(tmleOut_se)
methytmle <- sort(methytmle)
methytmle <- as.data.frame(methytmle)

################################################################################
### VOLCANO PLOT
################################################################################
tt_volcano <- methytmle %>%
  dplyr::arrange(pvalFDR) %>%
  dplyr::mutate(
    LATE = I(LocalATE),
    lowerLATE = ifelse(LocalATE > 0, lowerCI, upperCI),
    logPval = -log10(pvalue),
    cpgs = rownames(methytmle),
    top = c(rep(1, 50), rep(0, nrow(.) - 50)),
    colorMag = ifelse((LATE > 3.0) & (pvalFDR < 0.05), "1",
                       ifelse((LATE < -3.0) & (pvalFDR < 0.05), "-1", "0")),
    colorSign = ifelse(LATE > 0, "+1", "-1")
  ) %>%
  dplyr::select(which(colnames(.) %in% c("LowerLATE", "LATE", "logPval", "top",
                                         "colorMag", "colorSign", "cpgs"))) %>%
  na.omit()

pdf(file = here("graphs", "lateVolcPlot.pdf"))
p1 <- ggplot(tt_volcano, ggplot2::aes(x = LATE, y = logPval)) +
  geom_point(aes(colour = colorMag)) +
  geom_text(aes(label = ifelse(top != 0, as.character(cpgs), '')),
            hjust = 0, vjust = 0, check_overlap = TRUE) +
  xlab("Local Average Treatment Effect") +
  ylab("-log10(raw p-value)") +
  ggtitle("Local Average Treatment Effect vs. -log(p-value)") +
  xlim(c(min(tt_volcano$LATE) - 5, max(tt_volcano$LATE) + 15)) +
  scale_colour_manual(values = pal2[1:3], guide = FALSE) + theme_minimal()
print(p1)
dev.off()

pdf(file = here("graphs", "boundedVolcPlot.pdf"))
p2 <- ggplot(tt_volcano, ggplot2::aes(x = LowerATE, y = logPval)) +
  geom_point(aes(colour = colorSign)) +
  geom_text(aes(label = ifelse(top != 0, as.character(cpgs), '')),
            hjust = 0, vjust = 0, check_overlap = TRUE) +
  xlab("Lower Bound of Local Average Treatment Effect") +
  ylab("-log10(raw p-value)") +
  ggtitle("Local Average Treatment Effect (lower bound) vs. -log(p-value)") +
  xlim(c(min(tt_volcano$LowerATE), max(tt_volcano$LowerATE))) +
  scale_colour_manual(values = pal2[1:3], name = "L-ATE Direction",
                      labels = c("L-ATE > 0", "L-ATE < 0")) + theme_minimal()
print(p2)
dev.off()
################################################################################

################################################################################
### Manhattan plot
################################################################################
# first set up data
CHR <- as.factor(seqnames(subset(CpGs_gr,
                                 names(CpGs_gr) %in% rownames(methytmle))))
BP <- start(subset(CpGs_gr, names(CpGs_gr) %in% rownames(methytmle)))
inManhattan <- methytmle %>%
  dplyr::mutate(
    CHR = as.numeric(CHR),
    BP = as.numeric(BP),
    SNP = rownames(.),
    LATE = I(LocalATE),
    lowerLATE = ifelse(LocalATE > 0, lowerCI, upperCI),
    colorSign = ifelse(LATE > 0, "+1", "-1"),
    P = I(pvalue)
  ) %>%
  dplyr::select(which(colnames(.) %in% c("CHR", "BP", "SNP", "LATE", "P",
                                         "lowerLATE", "colorSign"))) %>%
  na.omit()

pdf(file = here("graphs", "CpG_ManhattPlot.pdf"))
manhattan(inManhattan,
          #p = "LATE",
          ylim = c(1, 15),
          cex = 0.6, cex.axis = 0.9,
          main = "Manhattan plot of subset of CpGs",
          #ylab = "LATE",
          col = c("blue4", "orange3"),
          suggestiveline = F,
          genomewideline = F,
          chrlabs = c(1:20, "P", "Q")
         )
dev.off()
################################################################################
