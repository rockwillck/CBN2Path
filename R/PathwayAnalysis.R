#' permutations
#'
#' @param n total number of elements in the set
#' @param r subset size
#' @param v 1:n
#' @param set Logical flag indicating whether duplicates should be removed from the source vector v. Defaults to TRUE.
#' @param repeatsAllowed Logical flag indicating whether the constructed vectors may include duplicated values. Defaults to FALSE.
#'
#' @return a matrix with (n!/(n-r)!) rows and r columns
#' @export
#'
#' @examples
#' perm <- permutations(4, 4)
permutations <- function(n, r, v = 1:n, set = TRUE, repeatsAllowed = FALSE) {
    if (mode(n) != "numeric" || length(n) != 1 || n < 1 || (n %% 1) != 0) {
        stop("bad value of n")
    }
    if (mode(r) != "numeric" || length(r) != 1 || r < 1 || (r %% 1) != 0) {
        stop("bad value of r")
    }
    if (!is.atomic(v) || length(v) < n) {
        stop("v is either non-atomic or too short")
    }
    if ((r > n) & repeatsAllowed == FALSE) {
        stop("r > n and repeats.allowed=FALSE")
    }
    if (set) {
        v <- unique(sort(v))
        if (length(v) < n) {
            stop("too few different elements")
        }
    }
    v0 <- vector(mode(v), 0)
    if (repeatsAllowed) {
        sub <- function(n, r, v) {
            if (r == 1) {
                matrix(v, n, 1)
            } else if (n == 1) {
                matrix(v, 1, r)
            } else {
                inner <- Recall(n, r - 1, v)
                cbind(rep(v, rep(nrow(inner), n)), matrix(t(inner),
                    ncol = ncol(inner), nrow = nrow(inner) * n,
                    byrow = TRUE
                ))
            }
        }
    } else {
        sub <- function(n, r, v) {
            if (r == 1) {
                matrix(v, n, 1)
            } else if (n == 1) {
                matrix(v, 1, r)
            } else {
                X <- NULL
                for (i in 1:n) {
                    X <- rbind(X, cbind(v[i], Recall(n -
                        1, r - 1, v[-i])))
                }
                X
            }
        }
    }
    sub(n, r, v[1:n])
}



#' generateMatrixGenotypes
#'
#' @param g genotype length
#'
#' @return a genotype matrix with ncol=g and nrow=2^g
#' @export
#'
#' @examples
#' generateMatrixGenotypes(4)
generateMatrixGenotypes <- function(g) {
    if (g > 20) {
        stop("This would generate more than one million genotypes")
    }
    f1 <- function(n) {
        lapply(seq.int(n), function(x) t(combn(n, x)))
    }
    genotNums <- f1(g)
    listOfVectors <- function(y) {
        lapply(unlist(lapply(y, function(x) {
            apply(x, 1, list)
        }), recursive = FALSE), function(m) m[[1]])
    }
    genotNums <- listOfVectors(genotNums)
    v <- rep(0, g)
    mat <- matrix(unlist(lapply(genotNums, function(x) {
        v[x] <- 1
        return(v)
    })), ncol = g, byrow = TRUE)
    mat <- rbind(rep(0, g), mat)
    colnames(mat) <- LETTERS[1:g]
    return(mat)
}



#' genotypeFeasibility
#'
#' @param genotypes the full set of potential binary genotypes of a given length.
#' @param dag matrix representing the DAG of restrictions.
#' @param x the number of mutations considered.
#'
#' @return a binary vector, which indicates feasibility or infeasibility of a set of genotypes
#' @export
#'
#' @examples
#' geno4 <- generateMatrixGenotypes(4)
#' dag <- matrix(c(4, 4, 4, 1, 2, 3), 3, 2)
#' x <- 4
#' genoF4 <- genotypeFeasibility(geno4, dag, x)
genotypeFeasibility <- function(genotypes, dag, x) {
    vec <- matrix(1, nrow = (2^x), ncol = 1)
    d <- dim(dag)[1]
    if (d > 0) {
        for (i in 1:(2^x)) {
            for (j in 1:d) {
                xx <- dag[j, 1]
                yy <- dag[j, 2]
                if ((genotypes[i, yy] == 1) && (genotypes[i, xx] == 0)) {
                    vec[i] <- 0
                }
            }
        }
    }
    return(vec)
}



#' pathwayFeasibility
#'
#' @param dag matrix representing the DAG of restrictions.
#' @param x the number of mutations considered.
#'
#' @return a binary vector, which indicates feasibility or infeasibility of a set of pathways
#' @export
#'
#' @examples
#' dag <- matrix(c(4, 4, 4, 1, 2, 3), 3, 2)
#' x <- 4
#' pathwayFeasibility(dag, x)
pathwayFeasibility <- function(dag, x) {
    perm <- permutations(x, x)
    p <- dim(perm)[1]
    d <- dim(dag)[1]
    vec <- numeric(p) + 1
    if (d > 0) {
        for (i in 1:p) {
            for (j in 1:d) {
                a1 <- dag[j, 1]
                a2 <- dag[j, 2]
                b1 <- which(perm[i, ] == a1)
                b2 <- which(perm[i, ] == a2)
                if (b2 < b1) {
                    vec[i] <- 0
                }
            }
        }
    }
    return(vec)
}




#' pathwayGenotypeCompatibility
#'
#' @param pathway a vector representing the given pathway.
#' @param genotype a binary vector representing the given genotype.
#'
#' @return returns 1 (if the given genotype is compatible with the given pathway), and 0 otherwise
#' @export
#'
#' @examples
#' geno1 <- c(1, 0, 1, 0)
#' geno2 <- c(1, 1, 0, 0)
#' path <- c(1, 2, 3, 4)
#' pathwayGenotypeCompatibility(path, geno1)
#' pathwayGenotypeCompatibility(path, geno2)
pathwayGenotypeCompatibility <- function(pathway, genotype) {
    c <- 1
    for (i in 1:(length(pathway) - 1)) {
        p <- pathway[i]
        if (genotype[p] == 0) {
            for (j in (i + 1):length(pathway)) {
                q <- pathway[j]
                if (genotype[q] == 1) {
                    c <- 0
                    break
                }
            }
        }
        if (c == 0) {
            break
        }
    }
    return(c)
}

#' base2IndVec
#'
#' @param vec a binary genotype vector
#'
#' @return a number used for indexing a given genotype
#' @export
#'
#' @examples
#' vec <- c(0, 1, 0, 1)
#' base2IndVec(vec)
base2IndVec <- function(vec) {
  tot <- 0
  for (i in 1:length(vec)) {
    if (vec[i] == 1) {
      tot <- tot + (2^i)
    }
  }
  return(tot)
}



#' pathwayCompatibilityQuartet
#'
#' @param gMat The n by 4 binary genotype matrix representing a given quartet for a sample of n genotypes.
#'
#' @return The compatibility score, which is represented as a vector of length 24, each element of which corresponds to one of the 24 pathways of length 4.
#' @export
#'
#' @examples
#' set.seed(100)
#' gMat <- matrix(sample(c(0, 1), 800, replace = TRUE), 200, 4)
#' pathwayCompatibilityQuartet(gMat)
pathwayCompatibilityQuartet <- function(gMat) {
  g <- generateMatrixGenotypes(4)
  indxs <- numeric(16)
  for (i in 1:16) {
    indxs[i] <- base2IndVec(g[i, ])
  }
  d <- dim(gMat)[1]
  v <- numeric(16)
  for (i in 1:d) {
    vec <- gMat[i, ]
    indx <- base2IndVec(vec)
    gIdx <- which(indxs == indx)
    v[gIdx] <- v[gIdx] + 1
  }
  v <- v / sum(v)
  perm <- permutations(4, 4)
  m <- matrix(0, 24, 16)
  for (i in 1:24) {
    for (j in 1:16) {
      m[i, j] <- pathwayGenotypeCompatibility(perm[i, ], g[j, ])
    }
  }
  c <- as.numeric(m %*% v)
  return(c)
}


#' genotypeMatrixMutator
#'
#' @param mat The genotype matrix including sampled genotypes, which need to be mutated.
#' @param fp False positive rate
#' @param fn False negative rate
#'
#' @return The mutated version of the genotype matrix
#' @export
#'
#' @examples
#' set.seed(100)
#' gMat <- matrix(sample(c(0, 1), 800, replace = TRUE), 200, 4)
#' gMatMut <- genotypeMatrixMutator(gMat, 0.2, 0.2)
genotypeMatrixMutator <- function(mat, fp, fn) {
  d <- dim(mat)[1]
  allP <- which(mat == 1)
  allN <- which(mat == 0)
  sampleFp <- allN[sample(1:length(allN), round(fp * length(allN)))]
  sampleFn <- allP[sample(1:length(allP), round(fn * length(allP)))]
  for (i in 1:length(sampleFp)) {
    x <- sampleFp[i] %% d
    y <- ceiling(sampleFp[i] / d)
    mat[x, y] <- 1
  }
  for (i in 1:length(sampleFn)) {
    x <- sampleFn[i] %% d
    y <- ceiling(sampleFn[i] / d)
    mat[x, y] <- 0
  }
  return(mat)
}


#' pathProbSSWM
#'
#' @param fitness A vector of length 2^x, each element of which representing the fitness assigned to one of the 2^x genotypes.
#' @param x The number of mutations considered.
#'
#' @return vector of probabilities assigned to all potential pathways of length x
#' @export
#'
#' @examples
#' f <- c(0, 0.1, 0.2, 0.1, 0.2, 0.4, 0.3, 0.2, 0.2, 0.1, 0, 0.6, 0.4, 0.3, 0.2, 1)
#' x <- 4
#' pathP <- pathProbSSWM(f, x)
pathProbSSWM <- function(fitness, x) {
  genotypes <- generateMatrixGenotypes(x)
  indx <- matrix(0, nrow = 2^x, ncol = 1)
  for (k in 1:(2^x)) {
    for (j in 1:x) {
      indx[k, 1] <- indx[k, 1] + 2^(j - 1) * genotypes[k, j]
    }
  }
  perm <- permutations(x, x)
  prob <- numeric(dim(perm)[1])
  tot <- 0
  for (i in 1:dim(perm)[1]) {
    temp1 <- 1
    vec <- perm[i, ]
    geno <- matrix(0, nrow = (x + 1), ncol = x)
    for (j in 1:x) {
      for (k in (j + 1):(x + 1)) {
        geno[k, (vec[j])] <- 1
      }
    }
    genoIndx <- matrix(0, nrow = (x + 1), ncol = 1)
    for (j in 1:(x + 1)) {
      for (k in 1:x) {
        genoIndx[j, 1] <- genoIndx[j, 1] + 2^(k - 1) * geno[j, k]
      }
    }
    fitnessVec <- matrix(0, nrow = (x + 1), ncol = 1)
    for (j in 1:(x + 1)) {
      fitnessVec[j] <- fitness[which(indx == genoIndx[j])]
    }
    flag <- 0
    for (j in 2:(x + 1)) {
      if (fitnessVec[j] < fitnessVec[(j - 1)]) {
        flag <- 1
      }
    }
    if (flag == 0) {
      for (j in 1:x) {
        sn <- which(geno[j, ] == 0)
        n <- length(sn)
        s <- fitnessVec[(j + 1)] - fitnessVec[j]
        t <- 0
        for (k in 1:n) {
          ggeno <- geno[j, ]
          ggeno[(sn[k])] <- 1
          ggenoIndx <- 0
          for (l in 1:x) {
            ggenoIndx <- ggenoIndx + 2^(l - 1) * ggeno[l]
          }
          fitness2 <- fitness[which(indx == ggenoIndx)]
          s1 <- fitness2 - fitnessVec[j]
          if (s1 > 0) {
            t <- t + s1
          }
        }
        temp1 <- temp1 * (s / t)
      }
      prob[i] <- temp1
    }
  }
  gg <- sum(prob, na.rm = TRUE)
  prob <- prob / gg
  prob <- as.numeric(prob)
  prob[which(is.na(prob))] <- 0
  return(prob)
}





#' pathProbCBN: quantifies pathway probabilities using the output of CT-CBN or H-CBN
#'
#' @param dag matrix representing the DAG of restrictions.
#' @param lambda the lambda values, which are produced by the CBN model.
#' @param x the number of mutations considered.
#'
#' @return vector of probabilities assigned to all potential pathways of length x
#' @export
#'
#' @examples
#' dag <- matrix(c(2, 2, 4, 1, 3, 3), 3, 2)
#' lambda <- c(1, 4, 3, 2.5, 2)
#' x <- 4
#' pathP <- pathProbCBN(dag, lambda, x)
pathProbCBN <- function(dag, lambda, x) {
    genotypes <- generateMatrixGenotypes(x)
    indx <- matrix(0, nrow = 2^x, ncol = 1)
    for (k in 1:(2^x)) {
        for (j in 1:x) {
            indx[k, 1] <- indx[k, 1] + 2^(j - 1) * genotypes[k, j]
        }
    }
    if (sum(is.na(dag)) > 0) { dag <- matrix(0, 0, 0) }
    allowedSet <- genotypeFeasibility(genotypes, dag, x)
    perm <- permutations(x, x)
    prob <- numeric(dim(perm)[1])
    tot <- 0
    for (i1 in 1:dim(perm)[1]) {
        temp1 <- 1
        vec <- perm[i1, ]
        geno <- matrix(0, nrow = (x + 1), ncol = x)
        for (j1 in 1:x) {
            for (k1 in (j1 + 1):(x + 1)) {
                geno[k1, (vec[j1])] <- 1
            }
        }
        flag <- 1
        for (j1 in 1:(x + 1)) {
            index <- 0
            for (k1 in 1:x) {
                index <- index + 2^(k1 - 1) * geno[j1, k1]
            }
            finalIndex <- which(indx == index)
            if (allowedSet[finalIndex] == 0) {
                flag <- 0
            }
        }
        for (j1 in 1:x) {
            set1 <- which(geno[j1, ] == 1)
            set2 <- which(geno[j1 + 1, ] == 1)
            indxLambda <- setdiff(set2, set1)
            sn <- which(geno[j1, ] == 0)
            sn2 <- numeric(length(sn))
            for (kaka in 1:length(sn)) {
                tempIndex <- 0
                for (kk in 1:x) {
                    tempIndex <- tempIndex + 2^(kk - 1) * geno[j1, kk]
                }
                tempIndex <- tempIndex + 2^(sn[kaka] - 1)
                finalIndex <- which(indx == tempIndex)
                sn2[kaka] <- allowedSet[finalIndex]
            }
            snn <- sn[which(sn2 == 1)]
            t <- sum(lambda[(snn + 1)])
            s <- lambda[(indxLambda + 1)]
            temp1 <- temp1 * (s / t)
        }
        if (flag == 0) {
            temp1 <- 0
        }
        prob[i1] <- temp1
    }
    prob <- as.numeric(prob)
    return(prob)
}


#' pathProbQuartetCTCBN
#'
#' @param gMat The n by 4 binary genotype matrix representing a given quartet for a sample of n genotypes.
#'
#' @return The probability distribution (returned by the CT-CBN model), which is represented as a vector of length 24.
#' @export
#'
#' @examples
#' set.seed(100)
#' gMat <- matrix(sample(c(0, 1), 12, replace = TRUE), 3, 4)
#' pathProbQuartetCTCBN(gMat)
pathProbQuartetCTCBN <- function(gMat) {
    posets <- readRDS(system.file("extdata", "Posets.rds", package = "CBN2Path"))
    bc <- Spock$new(
        poset = posets,
        numMutations = 4,
        genotypeMatrix = cbind(1, gMat)
    )
    results <- ctcbn(bc)
    logLik <- numeric(219)
    for (i in 1:219) {
        logLik[i] <- as.numeric(results[[i]]$summary[4])
    }
    indx <- which.max(logLik)
    lambda <- as.numeric(results[[indx]]$summary[5:9])
    if (indx == 1) {
        dag <- matrix(0, 0, 0)
    } else {
        dag <- posets[[indx]]
    }
    pathProb <- pathProbCBN(dag, lambda, 4)
    return(pathProb)
}



#' pathProbQuartetHCBN
#'
#' @param gMat The n by 4 binary genotype matrix representing a given quartet for a sample of n genotypes.
#'
#' @return The probability distribution (returned by the H-CBN model), which is represented as a vector of length 24.
#' @export
#'
#' @examples
#' set.seed(100)
#' gMat <- matrix(sample(c(0, 1), 12, replace = TRUE), 3, 4)
#' pathProbQuartetHCBN(gMat)
pathProbQuartetHCBN <- function(gMat) {
   posets <- readRDS(system.file("extdata", "Posets.rds", package = "CBN2Path"))
   bc <- Spock$new(
     poset = posets,
     numMutations = 4,
     genotypeMatrix = cbind(1, gMat)
   )
   results <- hcbn(bc)
   logLik <- numeric(219)
   for (i in 1:219) {
     logLik[i] <- as.numeric(results[[i]]$summary[4])
   }
   indx <- which.max(logLik)
   lambda <- as.numeric(results[[indx]]$summary[5:9])
   if (indx == 1) { dag <- matrix(0, 0, 0) }
   else { dag <- posets[[indx]] }
   pathProb <- pathProbCBN(dag, lambda, 4)
   return(pathProb)
}

#' posetWeightingRCBN
#'
#' @param vec The likelihood vector corresponding to a given set of posets
#'
#' @return The poset weight vector determined using the reciprocal ranking method
#' @export
#'
#' @examples
#' set.seed(100)
#' logLik <- runif(219)
#' w1 <- posetWeightingRCBN(logLik)
posetWeightingRCBN <- function(vec) {
    w <- numeric(length(vec))
    for (i in 1:length(vec)) {
        temp <- sort(vec, index.return = TRUE, decreasing = TRUE)$ix
        w[i] <- 1 / (which(temp == i))
    }
    w <- w / sum(w)
    return(w)
}

#' pathEdgeMapper
#'
#' @param x number of mutations to consider
#'
#' @return Pathway to edge compatibility matrix, each element of which indicates whether a given edge is included in the transitive closure of a given pathway (1) or not (0).
#' @export
#'
#' @examples
#' peMap <- pathEdgeMapper(4)
pathEdgeMapper <- function(x) {
    path <- permutations(x, x)
    edge <- permutations(x, 2)
    p <- dim(path)[1]
    e <- dim(edge)[1]
    peMap <- matrix(0, p, e)
    for (i in 1:p) {
        for (j in 1:e) {
            x1 <- which(path[i, ] == edge[j, 1])
            x2 <- which(path[i, ] == edge[j, 2])
            if (x1 < x2) {
                peMap[i, j] <- 1
            }
        }
    }
    return(peMap)
}


#' edgeMarginalized
#'
#' @param pathProb The pathway probabilities returned in the step 3 of the R-CBN algorithm
#' @param x        The number of mutations to consider
#'
#' @return returns the marginal probability of all the potential edges
#' @export
#'
#' @examples
#' dag <- matrix(c(2, 2, 4, 1, 3, 3), 3, 2)
#' lambda <- c(1, 4, 3, 2.5, 2)
#' x <- 4
#' pathP <- pathProbCBN(dag, lambda, x)
#' edgeMarginalized(pathP, x)
edgeMarginalized <- function(pathProb, x) {
    peMap <- pathEdgeMapper(x)
    d <- dim(peMap)[2]
    edgeProb <- numeric(d)
    for (i in 1:d) {
        indx <- which(peMap[, i] == 1)
        edgeProb[i] <- sum(pathProb[indx])
    }
    return(edgeProb)
}


#' pathwayWeightingRCBN
#'
#' @param edgeProb Marginal edge probabilities
#' @param peMap Pathway-edge compatibility matrix
#'
#' @return  The pathway weights (step 4 of the R-CBN algorithm)
#' @export
#'
#' @examples
#' dag <- matrix(c(2, 2, 4, 1, 3, 3), 3, 2)
#' lambda <- c(1, 4, 3, 2.5, 2)
#' x <- 4
#' pathP <- pathProbCBN(dag, lambda, x)
#' edgeProb <- edgeMarginalized(pathP, x)
#' peMap <- pathEdgeMapper(4)
#' pathwayWeightingRCBN(edgeProb, peMap)
pathwayWeightingRCBN <- function(edgeProb, peMap) {
    d <- dim(peMap)[1]
    w <- numeric(d)
    for (i in 1:d) {
        w[i] <- 1
        indx <- which(peMap[i, ] == 1)
        for (j in indx) {
            w[i] <- w[i] * edgeProb[j]
        }
    }
    w <- w / sum(w)
    return(w)
}



#' pathNormalization
#'
#' @param pathProb The pathway probabilities returned in the step 3 of the R-CBN algorithm
#' @param x        The number of mutations to consider
#'
#' @return  The updated pathway probabilities (the step 5 of the R-CBN algorithm)
#' @export
#'
#' @examples
#' dag <- matrix(c(2, 2, 4, 1, 3, 3), 3, 2)
#' lambda <- c(1, 4, 3, 2.5, 2)
#' x <- 4
#' pathP <- pathProbCBN(dag, lambda, x)
#' pathN <- pathNormalization(pathP, x)
pathNormalization <- function(pathProb, x) {
    peMap <- pathEdgeMapper(x)
    edgeProb <- edgeMarginalized(pathProb, x)
    w <- pathwayWeightingRCBN(edgeProb, peMap)
    pathProbn <- ((w * pathProb) / sum(w * pathProb))
    return(pathProbn)
}


#' pathProbQuartetRCBN
#'
#' @param gMat The n by 4 binary genotype matrix representing a given quartet for a sample of n genotypes.
#'
#' @return The probability distribution (returned by the R-CBN model), which is represented as a vector of length 24
#' @export
#'
#' @examples
#' set.seed(100)
#' gMat <- matrix(sample(c(0, 1), 12, replace = TRUE), 3, 4)
#' pathProbQuartetRCBN(gMat)
pathProbQuartetRCBN <- function(gMat) {
    posets <- readRDS(system.file("extdata", "Posets.rds", package = "CBN2Path"))
    bc <- Spock$new(
        poset = posets,
        numMutations = 4,
        genotypeMatrix = cbind(1, gMat)
    )
    results <- ctcbn(bc)
    logLik <- numeric(219)
    p <- matrix(0, 219, 24)
    for (i in 1:219) {
        logLik[i] <- as.numeric(results[[i]]$summary[4])
        lambda <- as.numeric(results[[i]]$summary[5:9])
        if (i == 1) {
            dag <- matrix(0, 0, 0)
        } else {
            dag <- posets[[i]]
        }
        p[i, ] <- pathProbCBN(dag, lambda, 4)
    }
    w1 <- posetWeightingRCBN(logLik)
    pathProb1 <- apply((w1 * p), 2, sum) / sum(w1)
    pathProb2 <- pathNormalization(pathProb1, 4)
    return(pathProb2)
}




#' base2Indexing
#'
#' @param mat A given poset represented by a binary matrix (in B-CBN)
#'
#' @return #Poset weight vectors based on the frequency of occurrence in the BCBN MCMC-sampling scheme.
#' @export
#'
#' @examples
#' set.seed(100)
#' mat <- matrix(sample(c(0, 1), 16, replace = TRUE), 4, 4)
#' base2Indexing(mat)
base2Indexing <- function(mat) {
    count <- 0
    num <- 0
    d <- dim(mat)[1]
    for (i in 1:d) {
        for (j in 1:d) {
            count <- count + 1
            if (mat[i, j] == 1) {
                num <- num + 2^count
            }
        }
    }
    return(num)
}


#' pathProbQuartetBCBN
#'
#' @param gMat The n by 4 binary genotype matrix representing a given quartet for a sample of n genotypes.
#'
#' @return The probability distribution (returned by the B-CBN model), which is represented as a vector of length 24.
#' @export
#'
#' @examples
#' set.seed(100)
#' gMat <- matrix(sample(c(0, 1), 12, replace = TRUE), 3, 4)
#' pathProbQuartetBCBN(gMat)
pathProbQuartetBCBN <- function(gMat) {
    genotypeMatrix <- cbind(1, gMat)
    posets <- readRDS(system.file("extdata", "Posets.rds", package = "CBN2Path"))
    bc <- Spock$new(
      poset = posets,
      numMutations = 4,
      genotypeMatrix = cbind(1, gMat)
    )
    results <- ctcbn(bc)
    p <- matrix(0, 219, 24)
    for (i in 1:219) {
        lambda <- as.numeric(results[[i]]$summary[5:9])
        if (i == 1) {
            dag <- matrix(0, 0, 0)
        } else {
            dag <- posets[[i]]
        }
        p[i, ] <- pathProbCBN(dag, lambda, 4)
    }
    posetSamples <- bcbn(gMat)
    posetIndex <- numeric(219)
    for (i in 2:219) {
        d <- dim(posets[[i]])[1]
        tempMat <- matrix(0, 4, 4)
        for (j in 1:d) {
            x1 <- posets[[i]][j, 1]
            x2 <- posets[[i]][j, 2]
            tempMat[x1, x2] <- 1
        }
        posetIndex[i] <- base2Indexing(tempMat)
    }
    wB <- numeric(219)
    for (i in 1:100000) {
        sampleIndex <- base2Indexing(posetSamples[[i]])
        indx <- which(posetIndex == sampleIndex)
        wB[indx] <- wB[indx] + 1
    }
    wB <- wB / sum(wB)
    pathProb <- apply((wB * p), 2, sum) / sum(wB)
    return(pathProb)
}





#' jensenShannonDivergence
#'
#' @param prob1 The first (discrete) probability distribution (vector)
#' @param prob2 The second (discrete) probability distribution (vector)
#'
#' @return Jensen Shannon Divergence between the two (discrete) probability distributions
#' @export
#'
#' @examples
#' set.seed(100)
#' gMat <- matrix(sample(c(0, 1), 12, replace = TRUE), 3, 4)
#' pathCT <- pathProbQuartetCTCBN(gMat)
#' pathH <- pathProbQuartetHCBN(gMat)
#' jensenShannonDivergence(pathCT, pathH)
jensenShannonDivergence <- function(prob1, prob2) {
    d <- 0
    for (i in 1:length(prob1)) {
        if (prob1[i] > 0) {
            d <- d + prob1[i] * log2(prob1[i] / (0.5 * prob1[i] + 0.5 * prob2[i]))
        }
        if (prob2[i] > 0) {
            d <- d + prob2[i] * log2(prob2[i] / (0.5 * prob1[i] + 0.5 * prob2[i]))
        }
    }
    dv <- (d / 2)
    return(dv)
}


#' predictability
#'
#' @param prob Pathway probability vector
#' @param x The length of genotype vectors
#'
#' @return predictability
#' @export
#'
#' @examples
#' set.seed(100)
#' gMat <- matrix(sample(c(0, 1), 12, replace = TRUE), 3, 4)
#' pathCT <- pathProbQuartetCTCBN(gMat)
#' pathH <- pathProbQuartetHCBN(gMat)
#' predC <- predictability(pathCT, 4)
#' predictability(pathH, 4)
predictability <- function(prob, x) {
    tot <- 0
    for (i in 1:length(prob)) {
        if (sum(prob[i], na.rm = TRUE) > 0) {
            tot <- tot - prob[i] * log(prob[i])
        }
    }
    pred <- 1 - (tot / log(factorial(x)))
    return(pred)
}

