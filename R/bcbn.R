#' B-CBN
#'
#' @param data Generated data
#' @param n_samples Number of samples <def: 25000>
#' @param theta Theta <def: 0>
#' @param epsilon Epsilon <def: 0.05>
#' @param n_chains N-Chains <def: 4>
#' @param thin Thin <def: 10>
#' @param n_cores Number of parallelized cores <def: 1>
#' @param Max_L The maximum number of iteration <def: 1000>
#'
#' @return A matrix
#' @export
#'
#' @examples
#' if (!require("rBCBN")) {
#'     install.packages(CBN2Path::getBCBNinstall(), repos = NULL, type = "source")
#' }
#' bcbn()
bcbn <- function(data = NULL, n_samples = 25000, theta = 0, epsilon = 0.05, n_chains = 4, thin = 10, Max_L =1000, n_cores = 1) {
  if (!requireNamespace("rBCBN", quietly = TRUE)) {
    stop("The bcbn function requires the rBCBN package. Install the package with:\n\n  install.packages(getBCBNinstall(), repos = NULL, type = \"source\")\n ")
  }
  if (!checkBCBNVersion()) {
    stop(sprintf("You need version %s of rBCBN. Install the correct version with:\n\n  install.packages(getBCBNinstall(), repos = NULL, type = \"source\")\n ", getBCBNVersion()))
  }
  if (is.null(data)) {
    rBCBN::bcbn_mcmc(n_samples = n_samples, theta = theta, epsilon = epsilon, n_chains = n_chains, thin = thin, Max_L = Max_L, n_cores = n_cores)
  } else {
    rBCBN::bcbn_mcmc(data, n_samples, theta, epsilon, n_chains, thin, Max_L, n_cores)
  }
}

#' Transitive Closure
#'
#' @param poset Poset matrix
#'
#' @return Poset matrix
#' @export
#'
#' @examples
#' poset <- matrix(0, 10, 10)
#'
#' poset[1, 2] <- 1
#' poset[2, 3] <- 1
#' poset[3, 4] <- 1
#' poset[5, 4] <- 1
#' poset[6, 7] <- 1
#' poset[8, 9] <- 1
#' poset[8, 10] <- 1
#' poset[6, 9] <- 1

#' transitive_closure(poset)
transitive_closure <- function(poset) {
  rBCBN::transitiveClosure(poset)
}

#' Generate Data
#'
#' @param poset Poset matrix
#' @param thetas Vector of theta values
#' @param eps Epsilon
#' @param N N
#'
#' @return A matrix
#' @export
#'
#' @examples
#' poset <- matrix(0, 10, 10)
#'
#' poset[1, 2] <- 1
#' poset[2, 3] <- 1
#' poset[3, 4] <- 1
#' poset[5, 4] <- 1
#' poset[6, 7] <- 1
#' poset[8, 9] <- 1
#' poset[8, 10] <- 1
#' poset[6, 9] <- 1
#'
#' tr <- transitive_closure(poset)
#' theta <- c(0.8, 0.7, 0.6, 0.7, 0.4, 0.25, 0.6, 0.75, 0.5, 0.2)
#' eps <- 0.1
#' N <- 400
#'
#' generate_data(tr, theta, eps, N)
generate_data <- function(poset, thetas, eps, N) {
  rBCBN::generateData(poset, thetas, eps, N)
}

