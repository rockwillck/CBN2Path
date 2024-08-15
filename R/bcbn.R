#' B-CBN
#'
#' @param data Generated data
#' @param n_samples Number of samples <def: 25000>
#' @param theta Theta <def: 0>
#' @param epsilon Epsilon <def: 0.05>
#' @param n_chains N-Chains <def: 4>
#' @param thin Thin <def: 10>
#' @param n_cores Number of parallelized cores <def: 1>
#'
#' @return A matrix
#' @export
#'
#' @examples
#' bcbn()
bcbn <-function(data = NULL, n_samples = 25000, theta = 0, epsilon = 0.05, n_chains = 4,thin = 10, n_cores = 1) {
  bcbn_mcmc(data, n_samples, theta, epsilon, n_chains, thin, n_cores)
}
