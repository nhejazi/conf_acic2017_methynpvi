# set up preliminaries to call the function
source(here("..", "lib", "optimPosit.R"))
source(here("..", "lib", "FDRmsa.R"))
source(here("..", "src", "2_methytmle_disc.R"))
source(here("..", "src", "3_methytmle_cont.R"))
library(foreach); library(parallel); library(doParallel)
library(SummarizedExperiment); library(tmle); library(tmle.npvi)
gLib = c("SL.mean", "SL.glm", "SL.bayesglm", "SL.knn", "SL.nnet")
QLib = c("SL.mean", "SL.glm", "SL.knn", "SL.gam", "SL.nnet")

# test the function based on regular TMLE
tmleOut <- methytmle_disc(sumExp = neuroIQse,
                          clusters = clusterCpGs,
                          outcomeVar = outcomesIQ,
                          targetSites = targetSites,
                          g_lib = gLib,
                          Q_lib = QLib
                         )

# test the function using TMLE-NPVI
npviOut <- methytmle_cont(sumExp = neuroIQse,
                          clusters = clusterCpGs,
                          outcomeVar = outcomesIQ,
                          targetSites = targetSites,
                          nullRange = c(-2, 2),
                          numberSites = 5
                         )
