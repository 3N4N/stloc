weightMatrix_nD = function(x, span = 0.5) {

  ncells = nrow(x)
  coords = as.matrix(x)
  d = as.matrix(dist(coords))

  # extract a weights vector per cell

  W_raw = sapply(seq_len(ncells), function(cell) {
    dvec = d[cell,]
    vals = rep(0, ncells)
    # vals[order(dvec)[1:ceiling(span*ncells)]] = seq(1, 0, length.out = ceiling(span*ncells))
    vals[order(dvec)[1:ceiling(span*ncells)]] = 1
    return(vals)
  }, simplify = FALSE)

  W = do.call(rbind, W_raw)

  return(W)

}

corTaylor <- function(x,w=1) {

  if(!inherits(x,"matrix")) {
    stop("Input must be inherit ’matrix’ class.")
  }

  if (length(w) == 1) {
    w <- rep(1, ncol(x))
  }

  x <- apply(x, 1, function(i) w*i)
  d <- ncol(x)

  return ((1/sqrt(d))*sd(eigen(cor(x),only.values=TRUE)$values))
}
