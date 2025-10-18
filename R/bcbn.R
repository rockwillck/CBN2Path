defaultData <- function() {
  poset <- matrix(0,3,3)
  poset[1,2] <- 1
  poset[2,3] <- 1
  tr <- transitiveClosure(poset)
  theta <- c(0.8, 0.7, 0.6)
  eps <- 0.1
  n <- 10
  generateData(tr, theta, eps, n)
}

#' B-CBN
#'
#' @param data Generated data
#' @param nSamples Number of samples <def: 25000>
#' @param theta Theta <def: 0>
#' @param epsilon Epsilon <def: 0.05>
#' @param nChains N-Chains <def: 4>
#' @param thin Thin <def: 10>
#' @param nCores Number of parallelized cores <def: 1>
#' @param maxL The maximum number of iteration <def: 1000>
#'
#' @return A matrix
#' @export
#'
#' @examples
#' bcbn()
bcbn <- function(data = defaultData(), nSamples = 25000, theta = 0, epsilon = 0.05, nChains = 4, thin = 10, maxL = 1000, nCores = 1) {
  if (nChains < nCores) {
    message(paste("Number of chains was less than number of cores. Using number of chains (", nChains, ") as thread count.", sep = ""))
  }
  n <- dim(data)[2]
  nCases <- dim(data)[1]
  mList <- list()
  edgeList <- list()
  l <- 0
  converged <- 0
  converged2 <- 0
  repeat {
    l <- l + 1

    retWorker <- function(i) {
      message(paste("chain:", i))
      message(theta)
      if (all(theta == 0)) { theta <- as.double(runif(n)) }
      if (length(mList) != 0) {
        edgesIn <- c(t(edgeList[[i]][nSamples][[1]]))
        theta <- as.double(mList[[i]][nSamples, 1:n])
        epsilon <- mList[[i]][nSamples, n + 1]
      } else {
        edgesIn <- as.integer(rep(0, n * n))
      }
      ret <- .C("sample_full_cbn_", theta, as.integer(n), as.double(epsilon), edgesIn, as.integer(nSamples), as.integer(thin), as.integer(c(t(data))), as.integer(nCases), thetaOut = as.double(rep(0, n * nSamples)), epsilonOut = as.double(rep(0, nSamples)), edgesOut = as.integer(rep(0, nSamples * n * n)), logPosteriorOut = as.double(rep(0, nSamples)))
    }

    if (exists("MulticoreParam", mode = "function") || exists("SnowParam", mode = "function")) {
      if(Sys.info()["sysname"] == "Windows") {
        p <- SnowParam(workers = min(nChains, nCores))
      } else {
        p <- MulticoreParam(workers = min(nChains, nCores))
      }
      rets <- bplapply(1:nChains, retWorker, BPPARAM = p)
    } else {
      message("Parallelization not found - running sequentially.")
      rets <- lapply(datasets, retWorker)
    }

    mList <- list()
    mcmcList <- list()
    edgeList <- list()
    for (i in 1:nChains) {
      thetaM <- matrix(rets[[i]]$thetaOut, ncol = n, byrow = TRUE)
      paraMatrix <- cbind(thetaM, rets[[i]]$epsilonOut, rets[[i]]$logPosteriorOut)
      mList[[i]] <- paraMatrix
      mcmcList[[i]] <- mcmc(paraMatrix)
      subList <- list()
      edgeSum <- 0
      for (k in 1:nSamples) {
        subList[[k]] <- matrix(rets[[i]]$edgesOut[((k - 1) * n * n + 1):(k * n * n)], ncol = n, byrow = TRUE)
        edgeSum <- edgeSum + matrix(rets[[i]]$edgesOut[((k - 1) * n * n + 1):(k * n * n)], ncol = n, byrow = TRUE)
      }
      edgeList[[i]] <- subList
      print(summary(paraMatrix))
    }
    mcList <- mcmc.list(mcmcList)
    terr <- try(gDiag <- gelman.diag(mcList))
    if (!inherits(terr, 'try-error')) {
      print(paste0("Criterion: ", max(gDiag$psrf[, 1])))
      print(gDiag)
      print("##########################################")
      if ((max(gDiag$psrf[, 1]) < 1.1) || (l > maxL)) {
        converged <- 1
        if (converged == 1) { break }
      }
      print(gDiag)
    }
    print(paste("finished run:", l))
  }
  eList <- c()
  for (i in 1:nChains) { eList <- c(eList, edgeList[[i]]) }
  return(eList)
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
#' poset[1, 2] <- 1
#' poset[2, 3] <- 1
#' poset[3, 4] <- 1
#' poset[5, 4] <- 1
#' poset[6, 7] <- 1
#' poset[8, 9] <- 1
#' poset[8, 10] <- 1
#' poset[6, 9] <- 1
#' transitiveClosure(poset)
transitiveClosure <- function(poset) {
  p <- nrow(poset)
  for (i in 2:p) {
    poset <- poset + matrixPower(poset, i)
  }
  poset <- apply(poset, 2, function(x) as.numeric(x > 0))
  return(poset)
}

#' Generate Data
#'
#' @param poset Poset matrix
#' @param theta Vector of theta values
#' @param eps Epsilon
#' @param n N
#'
#' @return A matrix
#' @export
#'
#' @examples
#' poset <- matrix(0, 10, 10)
#' poset[1, 2] <- 1
#' poset[2, 3] <- 1
#' poset[3, 4] <- 1
#' poset[5, 4] <- 1
#' poset[6, 7] <- 1
#' poset[8, 9] <- 1
#' poset[8, 10] <- 1
#' poset[6, 9] <- 1
#' tr <- transitiveClosure(poset)
#' theta <- c(0.8, 0.7, 0.6, 0.7, 0.4, 0.25, 0.6, 0.75, 0.5, 0.2)
#' eps <- 0.1
#' n <- 400
#' generateData(tr, theta, eps, n)
generateData <- function(poset, theta, eps, n) {
  nTheta <- length(theta)
  data <- matrix(0, n, nTheta)
  parentsList <- list()
  for (j in 1:nTheta) { parentsList[[j]] <- getParents(poset, j) }
  g <- graph.adjacency(poset)
  topoSorted <- topological.sort(g)
  for (i in 1:n) {
    for (j in topoSorted) {
      if (length(parentsList[[j]]) == 0 || prod(data[i, parentsList[[j]]]) == 1) {
        data[i, j] <- rbinom(1, 1, theta[j])
      }
    }
    epsNoise <- rbinom(nTheta, 1, eps)
    data[i, ] <- (data[i, ] + epsNoise) %% 2
  }
  data
}
