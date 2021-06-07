weightMatrix_nD = function(x, span = 0.5) {

  ncells = nrow(x)
  coords = as.matrix(x)
  d = as.matrix(dist(coords))

  # extract a weights vector per cell

  W_raw = sapply(seq_len(ncells), function(cell) {
    dvec = d[cell,]
    vals = rep(0, ncells)
    vals[order(dvec)[1:ceiling(span*ncells)]] = seq(1, 0, length.out = ceiling(span*ncells))
    # vals[order(dvec)[1:ceiling(span*ncells)]] = 1
    return(vals)
  }, simplify = FALSE)

  W = do.call(rbind, W_raw)

  return(W)

}

corTaylor <- function(x, w = 1) {

  if(!inherits(x, "matrix")) {
    stop("Input must be inherit ’matrix’ class.")
  }

  if (length(w) == 1) {
    w <- rep(1, ncol(x))
  }

  x <- apply(x, 1, function(i) w*i)
  d <- ncol(x)

  return ((1/sqrt(d)) * sd(eigen(cor(x),only.values=TRUE)$values))
}


maxEigenVal <- function(x, w=1) {

  if(!inherits(x,"matrix")) {
    stop("Input must be inherit ’matrix’ class.")
  }

  if (length(w) == 1) {
    w <- rep(1, ncol(x))
  }

  d <- nrow(x)
  x <- apply(x, 1, function(i) w*i)
  # cen_x <- apply(x, 1, function(i) i - ((w*i)/sum(w)))
  # cov_x <- (t(cen_x)%*%diag(w)%*%cen_x)/sum(w)
  # cor_x <- cov2cor(cov_x)

  # return (max(eigen(cor(x),only.values=TRUE)$values))
  # return (max(eigen(cov_x,only.values=TRUE)$values))
  # return (max(eigen(cor_x,only.values=TRUE)$values))
  return (max(eigen(t(x) %*% x,only.values=TRUE)$values))
}

zScore <- function(x, index, meanOfZenes, sdOfZenes) {

  zscore_for_single_gene = matrix(nrow = nrow(x), ncol = 1);
  sm = 0
  for(i in 1:(nrow(x))) {
    zscore_for_single_gene[i, 1] = (x[i, index] - meanOfZenes[i]) / sdOfZenes[i];
    sm = sm + zscore_for_single_gene[i, 1]
  }
  aggrZscore = sm / sqrt(nrow(x));
  return (aggrZscore)
}
