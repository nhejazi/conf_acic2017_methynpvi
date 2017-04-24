# select levels of A, W to discretize to while optimizing positivity assumption

optimPosit <- function(A, W, posMin = 0.1) {
  stopifnot(length(A) == nrow(W))

  nObs <- length(A)
  initGuessW <- round((posMin * nObs) / length(unique(A))) # heuristic W binning
  if (class(W) != "data.frame") W <- as.data.frame(W) # cover use of "ncol"
  outW <- NULL # concatenate W columnwise as we discretize each covar below

  for (obsW in 1:ncol(W)) {
    inW <- as.numeric(W[, obsW])
    discrW <- as.numeric(as.factor(gtools::quantcut(x = inW, q = initGuessW)))
    check <- sum((table(A, discrW) / nObs) < posMin)
    nextGuessW <- initGuessW
    while (check > 0) {
      nextGuessW <- (nextGuessW - 1)
      discrW <- as.numeric(as.factor(gtools::quantcut(x = inW, q = nextGuessW)))
      check <- sum((table(A, discrW) / nObs) < posMin)
    }
    outW <- cbind(outW, discrW)
  }
  out <- as.data.frame(outW)
  colnames(out) <- colnames(W)
  rownames(out) <- rownames(W)
  if(length(which(colSums(out) == nObs)) > 0) {
    out <- out[, -which(colSums(out) == nObs), drop = FALSE]
  }
  return(out)
}
