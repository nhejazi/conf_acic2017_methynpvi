library(here)
library(GenomicRanges)
library(parallel)
library(doParallel)
library(foreach)
library(limma)

makeFullData <- function() {

# setting project and data directories...
  proj_dir <- here()
  data_neuro <- here("data", "ch_450k_batch123_Neuro7yr_cr01_12162016OS.RDATA")
  data_m450k <- here("..", "..", "data", "holland-methylation",
                     "Batch123_ASMN_BMIQ_pp_ALL_04_NAkeep_XYkeep_mvalues_08142014.Rdata")

  # load data sets
  load(data_neuro)
  load(data_m450k)

  annot <- data_pp$annotation
  mvals <- data_pp$mmatrix

  # extract outcomes of interest
  outcomes <- c("fsiq_i_7y", "vciq_i_7y", "priq_i_7y", "wmiq_i_7y", "psiq_7y",
                "baattss_7y", "pbdss_7y", "ftdr_7y")

  # extract potential confounders
  confounders <- c("educcat", "sppstd_6m", "CD8T.bakulski", "CD4T.bakulski",
                   "NK.bakulski", "Bcell.bakulski", "Mono.bakulski",
                   "Gran.bakulski", "nRBC.bakulski")



}
