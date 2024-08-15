default_data <- function() {
  
  poset <- matrix(0,10,10)
  
  poset[1,2] <-1
  poset[2,3] <-1
  poset[3,4] <-1
  poset[5,4] <-1
  poset[6,7] <-1
  poset[8,9] <-1
  poset[8,10] <-1
  poset[6,9] <-1
  
  tr<-transitiveClosure(poset)
  theta <- c(0.8, 0.7, 0.6, 0.7, 0.4, 0.25, 0.6, 0.75, 0.5, 0.2)
  eps <- 0.1
  N <- 400
  
  generateData(tr, theta, eps, N)
}

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
bcbn <-function(data = default_data(), n_samples = 25000, theta = 0, epsilon = 0.05, n_chains = 4,thin = 10, n_cores = 1) {
  bcbn_mcmc(data, n_samples, theta, epsilon, n_chains, thin, n_cores)
}
