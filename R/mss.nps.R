#' Minimum sample size using ALC with HPD interval for estimating the NPS
#'
#' @param len.max maximum admissible length for the HPD interval.
#' @param rho A number in (0, 1). The probability associated to the HPD interval is \code{1-rho}.
#' @param a vector parameter of the prior Dirichlet distribution.
#' @param L Number of replicates used in the simulation to estimate the HPD interval length. Default is 1000.
#' @param N Number of replicates used in the simulation of the posterior distribution of the NPS. Default is 100.
#' @param n0 A positive integer representing the initial sample size in which the function will check the criterion. Default is 1.
#'
#' @return Minimum sample size for estimating the NPS via the HPD interval.
#' @export
#'
#' @examples mss.nps(len.max = 0.5)
mss.nps <- function(len.max, rho = 0.05, a = rep(1, 3), L = 1E3, N = 1E2, n0 = 1) {
  packageStartupMessage("Computing...")
  cl <- match.call()
  n <- n0
  len <- len.max + 1
  while (mean(len) > len.max) {
    n <- n + 1
    len <- numeric(L)
    for (i in 1:L) {
      theta <- rdiri(N = 1, alpha = a)
      x <- as.vector(rmultinom(n = 1, size = n, prob = theta))
      NPS.pos <- rpos.nps(n = N, x = x, a = a)
      ab <- TeachingDemos::emp.hpd(NPS.pos, conf = 1 - rho)
      len[i] <- ab[2] - ab[1]
    }
  }
  # Output
  cat("\nCall:\n")
  print(cl)
  cat("\nMinimum Sample Size for NPS:\n")
  cat("n  = ", n, "\n")
}
