#' Random generation for the Dirichlet distribution
#'
#' @param N number of samples to generate.
#' @param alpha vector parameter.
#'
#' @return A matrix with the generated samples.
#'
#'
#' @examples rdiri(N = 100, alpha = rep(1, 3))
rdiri <- function(N = 1, alpha){
  out <- matrix(rgamma(N*length(alpha), shape = alpha), nrow = N,
                ncol = length(alpha), byrow = TRUE)
  out <- out/apply(out, 1, sum)
  return(out)
}
