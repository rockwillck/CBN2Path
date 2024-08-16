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
#' \dontrun{
#' bcbn()
#' }
bcbn <-function(data = NULL, n_samples = 25000, theta = 0, epsilon = 0.05, n_chains = 4,thin = 10, n_cores = 1) {
  if (is.null(data)){
    bcbn_mcmc(n_samples = n_samples, theta = theta, epsilon = epsilon, n_chains = n_chains, thin = thin, n_cores = n_cores)
  } else {
    bcbn_mcmc(data, n_samples, theta, epsilon, n_chains, thin, n_cores)
  }
}
