#' Random generation for the posterior distribution of the NPS
#'
#' @param n number of observations for random generation.
#' @param x a vector with the observed counts of detractors, passives and promoters, respectively.
#' @param a vector parameter of the prior Dirichlet distribution.
#'
#' @return A vector with a random sample from the posterior distribution of the NPS.
#' @export
#'
#' @examples rpos.nps(n = 100, x = c(10, 50, 80))
rpos.nps <- function(n, x, a = rep(1, 3)) {
  cl <- match.call()
  if (is.matrix(x)) {
    s <- colSums(x)
  }
  if (is.vector(x)) {
    s <- x
  }
  a.s <- a + s
  theta.pos <- rdiri(N = n, alpha = a.s)
  NPS.pos <- theta.pos[, 3] - theta.pos[, 1]
  return(NPS.pos)
}
