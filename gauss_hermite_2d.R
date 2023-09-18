library(statmod)

## find nodes and weights of a 2d Gaussian-Hermite quadrature
gauss_hermite_2d <- function(nn, mu, L, prune = 0.3){
  n1 <- nn[1]
  n2 <- nn[2]
  pts1 <- as.data.frame(statmod::gauss.quad(n1, kind="hermite"))
  pts2 <- as.data.frame(statmod::gauss.quad(n2, kind="hermite"))
  nodes0 <- expand.grid(pts1$nodes, pts2$nodes)
  nodes <- t( mu + sqrt(2) * L %*% t(nodes0) )
  w <- expand.grid(pts1$weights, pts2$weights)
  weights <- apply(w, 1, prod) * (1/pi)
  ## prune
  if(!is.null(prune)) {
    qwt <- quantile(weights, probs=prune)
    nodes <- nodes[weights > qwt,]
    weights <- weights[weights > qwt]
  }
  return(list(nodes = nodes, weights = weights))
}
