#' Bayesian inference for the NPS
#'
#' @param x a vector with the observed counts of detractors, passives and promoters, respectively.
#' @param a vector parameter of the prior Dirichlet distribution.
#' @param rho A number in (0, 1). The probability associated to the HPD interval is \code{1-rho}.
#' @param N Number of replicates used in the simulation. Default is 1000.
#'
#' @return Bayesian summaries for the posterior distribution of the NPS as mean, standard deviation and HPD interval.
#' @export
#'
#' @examples bayes.nps(x = c(10, 50, 80))
bayes.nps <- function(x, a = rep(1, 3), rho = 0.05, N = 1E3) {
  cl <- match.call()
  NPS.pos <- rpos.nps(n = N, x = x, a = a)
  ab <- TeachingDemos::emp.hpd(NPS.pos, conf = 1 - rho)
  # Output
  cat("\nCall:\n")
  print(cl)
  cat("\nBayesian Inference for NPS:\n")
  print(matrix(c(mean(NPS.pos), sd(NPS.pos), ab[1], ab[2]), 1, 4,
               dimnames = list("NPS", c("mean", "sd", "l-HPD", "u-HPD"))))
  cat("\n")
}
