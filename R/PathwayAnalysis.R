##############################################################################################################
#' combinations
#'
#' @param n total number of elements in the set
#' @param r subset size
#' @param v 1:n
#' @param set
#' @param repeats.allowed
#'
#' @return a matrix with (n choose r) rows and r columns
#' @export
#'
#' @examples
#' COMB<-combinations(10,4)
combinations <- function(n, r, v = 1:n, set = TRUE, repeats.allowed = FALSE) {
    if (mode(n) != "numeric" || length(n) != 1 || n < 1 || (n %% 1) !=
        0) {
        stop("bad value of n")
    }
    if (mode(r) != "numeric" || length(r) != 1 || r < 1 || (r %% 1) !=
        0) {
        stop("bad value of r")
    }
    if (!is.atomic(v) || length(v) < n) {
        stop("v is either non-atomic or too short")
    }
    if ((r > n) & repeats.allowed == FALSE) {
        stop("r > n and repeats.allowed=FALSE")
    }
    if (set) {
        v <- unique(sort(v))
        if (length(v) < n) {
            stop("too few different elements")
        }
    }
    v0 <- vector(mode(v), 0)
    if (repeats.allowed) {
        sub <- function(n, r, v) {
            if (r == 0) {
                v0
            } else if (r == 1) {
                matrix(v, n, 1)
            } else if (n == 1) {
                matrix(v, 1, r)
            } else {
                rbind(cbind(v[1], Recall(n, r - 1, v)), Recall(n -
                    1, r, v[-1]))
            }
        }
    } else {
        sub <- function(n, r, v) {
            if (r == 0) {
                v0
            } else if (r == 1) {
                matrix(v, n, 1)
            } else if (r == n) {
                matrix(v, 1, n)
            } else {
                rbind(
                    cbind(v[1], Recall(n - 1, r - 1, v[-1])),
                    Recall(n - 1, r, v[-1])
                )
            }
        }
    }
    sub(n, r, v[1:n])
}



#' permutations
#'
#' @param n total number of elements in the set
#' @param r subset size
#' @param v 1:n
#' @param set
#' @param repeats.allowed
#'
#' @return
#' @export
#'
#' @examples
#' PERM<-permutations(4,4)
permutations <- function(n, r, v = 1:n, set = TRUE, repeats.allowed = FALSE) {
    if (mode(n) != "numeric" || length(n) != 1 || n < 1 || (n %% 1) !=
        0) {
        stop("bad value of n")
    }
    if (mode(r) != "numeric" || length(r) != 1 || r < 1 || (r %% 1) !=
        0) {
        stop("bad value of r")
    }
    if (!is.atomic(v) || length(v) < n) {
        stop("v is either non-atomic or too short")
    }
    if ((r > n) & repeats.allowed == FALSE) {
        stop("r > n and repeats.allowed=FALSE")
    }
    if (set) {
        v <- unique(sort(v))
        if (length(v) < n) {
            stop("too few different elements")
        }
    }
    v0 <- vector(mode(v), 0)
    if (repeats.allowed) {
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



#' generate_matrix_genotypes
#'
#' @param g genotype length
#'
#' @return a genotype matrix with ncol=g and nrow=2^g
#' @export
#'
#' @examples
#' Geno4<-generate_matrix_genotypes(4)
generate_matrix_genotypes <- function(g) {
    if (g > 20) {
        stop("This would generate more than one million genotypes")
    }
    f1 <- function(n) {
        lapply(seq.int(n), function(x) combinations(n = n, r = x))
    }
    genotNums <- f1(g)
    list.of.vectors <- function(y) {
        lapply(unlist(lapply(y, function(x) {
            apply(x, 1, list)
        }), recursive = FALSE), function(m) m[[1]])
    }
    genotNums <- list.of.vectors(genotNums)
    v <- rep(0, g)
    mat <- matrix(unlist(lapply(genotNums, function(x) {
        v[x] <- 1
        return(v)
    })), ncol = g, byrow = TRUE)
    mat <- rbind(rep(0, g), mat)
    colnames(mat) <- LETTERS[1:g]
    return(mat)
}



#' Genotype_Feasibility
#'
#' @param genotypes the full set of potential binary genotypes of a given length.
#' @param DAG matrix representing the DAG of restrictions.
#' @param x the number of mutations considered.
#'
#' @return a binary vector, which indicates feasibility or infeasibility of a set of genotypes
#' @export
#'
#' @examples
#' Geno4<-generate_matrix_genotypes(4)
#' DAG<-matrix(c(4,4,4,1,2,3),3,2)
#' x<-4
#' GenoF4<-Genotype_Feasibility(Geno4,DAG,x)
Genotype_Feasibility <- function(genotypes, DAG, x) {
    vec <- matrix(1, nrow = (2^x), ncol = 1)
    D <- dim(DAG)[1]
    if (D > 0) {
        for (i in 1:(2^x)) {
            for (j in 1:D) {
                xx <- DAG[j, 1]
                yy <- DAG[j, 2]
                if ((genotypes[i, yy] == 1) && (genotypes[i, xx] == 0)) {
                    vec[i] <- 0
                } ## if the mutation ordering is not respected the genotype is labelled as infeasible.
            }
        }
    }
    return(vec)
}



#' Pathway_Feasibility
#'
#' @param DAG matrix representing the DAG of restrictions.
#' @param x the number of mutations considered.
#'
#' @return a binary vector, which indicates feasibility or infeasibility of a set of pathways
#' @export
#'
#' @examples
#' DAG<-matrix(c(4,4,4,1,2,3),3,2)
#' x<-4
#' PathF<-Pathway_Feasibility(DAG, x)
Pathway_Feasibility <- function(DAG, x) {
    PERM <- permutations(x, x) ### all x! possible permutations (mutational pathways)
    P <- dim(PERM)[1]
    D <- dim(DAG)[1]
    vec <- numeric(P) + 1
    if (D > 0) {
        for (i in 1:P) {
            for (j in 1:D) {
                a1 <- DAG[j, 1]
                a2 <- DAG[j, 2]
                b1 <- which(PERM[i, ] == a1)
                b2 <- which(PERM[i, ] == a2)
                if (b2 < b1) {
                    vec[i] <- 0
                }
            }
        }
    }
    return(vec)
}




#' Pathway_Genotype_Compatiblility
#'
#' @param Pathway a vector representing the given pathway.
#' @param Genotype a binary vector representing the given genotype.
#'
#' @return returns 1 (if the given genotype is compatible with the given pathway), and 0 otherwise
#' @export
#'
#' @examples
#' Geno1<-c(1,0,1,0)
#' Geno2<-c(1,1,0,0)
#' Path<-c(1,2,3,4)
#' Pathway_Genotype_Compatiblility(Path,Geno1)
#' Pathway_Genotype_Compatiblility(Path,Geno2)
Pathway_Genotype_Compatiblility <- function(Pathway, Genotype) {
    C <- 1 # by default compatible unless:
    for (i in 1:(length(Pathway) - 1)) {
        P <- Pathway[i]
        if (Genotype[P] == 0) {
            for (j in (i + 1):length(Pathway)) {
                Q <- Pathway[j]
                if (Genotype[Q] == 1) {
                    C <- 0
                    break
                }
            }
        }
        if (C == 0) {
            break
        }
    }
    return(C)
}



#' PathProb_CBN: quantifies pathway probabilities using the output of CT-CBN or H-CBN
#'
#' @param DAG matrix representing the DAG of restrictions.
#' @param LAMBDA the lambda values, which are produced by the CBN model.
#' @param x the number of mutations considered.
#'
#' @return vector of probabilities assigned to all potential pathways of length x
#' @export
#'
#' @examples
#' DAG<-matrix(c(2,2,4,1,3,3),3,2)
#' LAMBDA<-c(1,4,3,2.5,2)
#' x<-4
#' PathP<-PathProb_CBN(DAG, LAMBDA, x)
PathProb_CBN <- function(DAG, LAMBDA, x) {
    ### Step1: genotypes
    genotypes <- generate_matrix_genotypes(x) ## generates the genotype space[requires "OncoSimulR" package--two lines above]
    indx <- matrix(0, nrow = 2^x, ncol = 1) ## indexing the genotypes for easier retrival
    for (k in 1:(2^x)) {
        for (j in 1:x) {
            indx[k, 1] <- indx[k, 1] + 2^(j - 1) * genotypes[k, j]
        }
    }

    ### Step2: Allowed genotypes according to the DAG of restrictions (DAG)
    allowed_set <- Genotype_Feasibility(genotypes, DAG, x)

    ### Step3: Pathway Probabilities
    PERM <- permutations(x, x) ## all x! possible permutations (mutational pathways)
    Prob <- numeric(dim(PERM)[1]) ## pathway probabilities
    TOT <- 0 ## (the normalization factor)
    for (i1 in 1:dim(PERM)[1]) {
        TEMP1 <- 1 # temporarily stores the probability (later needs to be normalized)
        vec <- PERM[i1, ] # (i1)-th pathway
        GENO <- matrix(0, nrow = (x + 1), ncol = x) # constructing the possible x+1 genotypes, each corresponding to a given step of the given mutational pathway.
        for (j1 in 1:x) {
            for (k1 in (j1 + 1):(x + 1)) {
                GENO[k1, (vec[j1])] <- 1
            }
        }

        flag <- 1
        for (j1 in 1:(x + 1)) {
            index <- 0 # index of the (j1-th) genotype in the mutational pathway.
            for (k1 in 1:x) {
                index <- index + 2^(k1 - 1) * GENO[j1, k1]
            }
            final_index <- which(indx == index) # finding the index of the genotype in the genotype space
            if (allowed_set[final_index] == 0) {
                flag <- 0
            } # to check whether each genotype visited is allowed.
        }

        for (j1 in 1:x) {
            set1 <- which(GENO[j1, ] == 1)
            set2 <- which(GENO[j1 + 1, ] == 1)
            indx_lambda <- setdiff(set2, set1) # the index of the (j1-th) genotype in the pathway to retrieve the corresponding Lambda value
            SN <- which(GENO[j1, ] == 0) # the set of genotypes with one additional mutation than the current genotype
            ###################################################
            SN2 <- numeric(length(SN))
            for (kaka in 1:length(SN)) {
                TEMP_index <- 0
                for (kk in 1:x) {
                    TEMP_index <- TEMP_index + 2^(kk - 1) * GENO[j1, kk]
                }
                TEMP_index <- TEMP_index + 2^(SN[kaka] - 1)
                FINAL_index <- which(indx == TEMP_index)

                SN2[kaka] <- allowed_set[FINAL_index]
            }
            SNN <- SN[which(SN2 == 1)] # THE EXIT SET: the set of (allowed) genotypes with one additional mutation than the current genotype
            T <- sum(LAMBDA[(SNN + 1)]) # sum of the lambdas of the exit set [The denominator of the equation (10) in the main text]
            ###################################################
            S <- LAMBDA[(indx_lambda + 1)] # the lambda of the (j1-th mutation) [The numerator of the equation (10) in the main text]
            TEMP1 <- TEMP1 * (S / T) # The multiplication in the equation (10) in the main text
        }
        if (flag == 0) {
            TEMP1 <- 0
        } # If the pathway is infeasible, its probability will be zero.
        Prob[i1] <- TEMP1 # pathway probability
    }
    Prob <- as.numeric(Prob)
    return(Prob)
}


#' PathProb_Quartet_CTCBN
#'
#' @param gMat The n by 4 binary genotype matrix representing a given quartet for a sample of n genotypes.
#'
#' @return The probability distribution (returned by the CT-CBN model), which is represented as a vector of length 24.
#' @export
#'
#' @examples
#' set.seed(100)
#' gMat<-matrix(sample(c(0,1),800,replace = TRUE),200,4)
#' PathCT<-PathProb_Quartet_CTCBN(gMat)
PathProb_Quartet_CTCBN <- function(gMat) {
    Posets <- readRDS(system.file("extdata", "Posets.rds", package = "CBN2Path"))
    bc <- Spock$new(
        poset = Posets,
        numMutations = 4,
        genotypeMatrix = cbind(1, gMat)
    )
    Results <- ctcbn(bc)
    LogLik <- numeric(219)
    for (i in 1:219) {
        LogLik[i] <- as.numeric(Results[[i]]$summary[4])
    }
    INDX <- which.max(LogLik)
    LAMBDA <- as.numeric(Results[[INDX]]$summary[5:9])

    if (INDX == 1) {
        DAG <- matrix(0, 0, 0)
    } else {
        DAG <- Posets[[INDX]]
    }
    PathProb <- PathProb_CBN(DAG, LAMBDA, 4)
    return(PathProb)
}



#' PathProb_Quartet_HCBN
#'
#' @param gMat The n by 4 binary genotype matrix representing a given quartet for a sample of n genotypes.
#'
#' @return The probability distribution (returned by the H-CBN model), which is represented as a vector of length 24.
#' @export
#'
#' @examples
#' set.seed(100)
#' gMat<-matrix(sample(c(0,1),800,replace = TRUE),200,4)
#' PathH<-PathProb_Quartet_HCBN(gMat)
PathProb_Quartet_HCBN<-function(gMat){
   Posets<-readRDS(system.file("extdata","Posets.rds",package="CBN2Path"))
   bc=Spock$new(
     poset=Posets,
     numMutations=4,
     genotypeMatrix=cbind(1,gMat)
   )
   Results<-hcbn(bc)
   LogLik<-numeric(219)
   for (i in 1:219){
     LogLik[i]<-as.numeric(Results[[i]]$summary[4])
   }
   INDX<-which.max(LogLik)
   LAMBDA<-as.numeric(Results[[INDX]]$summary[5:9])

   if (INDX==1){DAG<-matrix(0,0,0)}
   else {DAG<-Posets[[INDX]]}
   PathProb<-PathProb_CBN(DAG,LAMBDA,4)
   return(PathProb)
}

#' Poset_Weighting_RCBN
#'
#' @param vec The likelihood vector corresponding to a given set of posets
#'
#' @return The poset weight vector determined using the reciprocal ranking method
#' @export
#'
#' @examples
#' set.seed(100)
#' LogLik<-runif(219)
#' W1<-Poset_Weighting_RCBN(LogLik)
Poset_Weighting_RCBN <- function(vec) {
    w <- numeric(length(vec))
    for (i in 1:length(vec)) {
        temp <- sort(vec, index.return = TRUE, decreasing = TRUE)$ix
        w[i] <- 1 / (which(temp == i))
    }
    w <- w / sum(w)
    return(w)
}

#' Path_Edge_Mapper
#'
#' @param x number of mutations to consider
#'
#' @return Pathway to edge compatibility matrix, each element of which indicates whether a given edge is included in the transitive closure of a given pathway (1) or not (0).
#' @export
#'
#' @examples
#' PEmap<-Path_Edge_Mapper(4)
Path_Edge_Mapper <- function(x) {
    PATH <- permutations(x, x)
    EDGE <- permutations(x, 2)
    P <- dim(PATH)[1]
    E <- dim(EDGE)[1]
    PEmap <- matrix(0, P, E)
    for (i in 1:P) {
        for (j in 1:E) {
            x1 <- which(PATH[i, ] == EDGE[j, 1])
            x2 <- which(PATH[i, ] == EDGE[j, 2])
            if (x1 < x2) {
                PEmap[i, j] <- 1
            }
        }
    }
    return(PEmap)
}


#' EdgeMarginalized
#'
#' @param PathProb The pathway probabilities returned in the step 3 of the R-CBN algorithm
#' @param x
#'
#' @return returns the marginal probability of all the potential edges
#' @export
#'
#' @examples
#' DAG<-matrix(c(2,2,4,1,3,3),3,2)
#' LAMBDA<-c(1,4,3,2.5,2)
#' x<-4
#' PathP<-PathProb_CBN(DAG, LAMBDA, x)
#' EdgeProb<-EdgeMarginalized(PathP,x)
EdgeMarginalized <- function(PathProb, x) {
    PEmap <- Path_Edge_Mapper(x)
    D <- dim(PEmap)[2]
    EdgeProb <- numeric(D)
    for (i in 1:D) {
        INDX <- which(PEmap[, i] == 1)
        EdgeProb[i] <- sum(PathProb[INDX])
    }
    return(EdgeProb)
}


#' Pathway_Weighting_RCBN
#'
#' @param EdgeProb Marginal edge probabilities
#' @param PEmap Pathway-edge compatibility matrix
#'
#' @return  The pathway weights (step 4 of the R-CBN algorithm)
#' @export
#'
#' @examples
#' DAG<-matrix(c(2,2,4,1,3,3),3,2)
#' LAMBDA<-c(1,4,3,2.5,2)
#' x<-4
#' PathP<-PathProb_CBN(DAG, LAMBDA, x)
#' EdgeProb<-EdgeMarginalized(PathP,x)
#' PEmap<-Path_Edge_Mapper(4)
#' W2<-Pathway_Weighting_RCBN(EdgeProb,PEmap)
Pathway_Weighting_RCBN <- function(EdgeProb, PEmap) {
    D <- dim(PEmap)[1]
    w <- numeric(D)
    for (i in 1:D) {
        w[i] <- 1
        INDX <- which(PEmap[i, ] == 1)
        for (j in INDX) {
            w[i] <- w[i] * EdgeProb[j]
        }
    }
    w <- w / sum(w)
    return(w)
}



#' Path_Normalization
#'
#' @param PathProb The pathway probabilities returned in the step 3 of the R-CBN algorithm
#' @param x
#'
#' @return  The updated pathway probabilities (the step 5 of the R-CBN algorithm)
#' @export
#'
#' @examples
#' DAG<-matrix(c(2,2,4,1,3,3),3,2)
#' LAMBDA<-c(1,4,3,2.5,2)
#' x<-4
#' PathP<-PathProb_CBN(DAG, LAMBDA, x)
#' PathN<-Path_Normalization(PathP, x)
Path_Normalization <- function(PathProb, x) {
    ### Step 4 of the R-CBN algorithm
    PEmap <- Path_Edge_Mapper(x)
    EdgeProb <- EdgeMarginalized(PathProb, x)
    w <- Pathway_Weighting_RCBN(EdgeProb, PEmap)
    ### Step 5 of the R-CBN algorithm
    PathProbn <- ((w * PathProb) / sum(w * PathProb)) # The normalized pathway probability
    return(PathProbn)
}


#' PathProb_Quartet_RCBN
#'
#' @param gMat The n by 4 binary genotype matrix representing a given quartet for a sample of n genotypes.
#'
#' @return The probability distribution (returned by the R-CBN model), which is represented as a vector of length 24
#' @export
#'
#' @examples
#' set.seed(100)
#' gMat<-matrix(sample(c(0,1),800,replace = TRUE),200,4)
#' PathR<-PathProb_Quartet_RCBN(gMat)
PathProb_Quartet_RCBN <- function(gMat) {
    ### Step 1: Constructing the P matrix
    Posets <- readRDS(system.file("extdata", "Posets.rds", package = "CBN2Path"))
    bc <- Spock$new(
        poset = Posets,
        numMutations = 4,
        genotypeMatrix = cbind(1, gMat)
    )
    Results <- ctcbn(bc)
    LogLik <- numeric(219)
    P <- matrix(0, 219, 24)
    for (i in 1:219) {
        LogLik[i] <- as.numeric(Results[[i]]$summary[4])
        LAMBDA <- as.numeric(Results[[i]]$summary[5:9])
        if (i == 1) {
            DAG <- matrix(0, 0, 0)
        } else {
            DAG <- Posets[[i]]
        }
        P[i, ] <- PathProb_CBN(DAG, LAMBDA, 4)
    }
    ### Step 2: Poset-Level weighting
    w1 <- Poset_Weighting_RCBN(LogLik)
    ### Step 3: Aggregating the Probability Distributions
    PathProb1 <- apply((w1 * P), 2, sum) / sum(w1)
    ### Step 4: Pathway-level weighting && Step 5: Updating the pathway probabilities
    PathProb2 <- Path_Normalization(PathProb1, 4)

    return(PathProb2)
}




#' Base2Indexing
#'
#' @param mat A given poset represented by a binary matrix (in B-CBN)
#'
#' @return #Poset weight vectors based on the frequency of occurence in the BCBN MCMC-sampling scheme.
#' @export
#'
#' @examples
#' set.seed(100)
#' mat<-matrix(sample(c(0,1),16,replace=TRUE),4,4)
#' Index<-Base2Indexing(mat)
Base2Indexing <- function(mat) {
    count <- 0
    num <- 0
    D <- dim(mat)[1]
    for (i in 1:D) {
        for (j in 1:D) {
            count <- count + 1
            if (mat[i, j] == 1) {
                num <- num + 2^count
            }
        }
    }
    return(num)
}


#' PathProb_Quartet_BCBN
#'
#' @param gMat The n by 4 binary genotype matrix representing a given quartet for a sample of n genotypes.
#'
#' @return The probability distribution (returned by the B-CBN model), which is represented as a vector of length 24.
#' @export
#'
#' @examples
#' set.seed(100)
#' gMat<-matrix(sample(c(0,1),800,replace = TRUE),200,4)
#' PathB<-PathProb_Quartet_BCBN(gMat)
PathProb_Quartet_BCBN <- function(gMat) {
    ### Step 1: Constructing the P matrix
    genotypeMatrix <- cbind(1, gMat)
    Posets <- readRDS(system.file("extdata", "Posets.rds", package = "CBN2Path"))
    bc <- Spock$new(
      poset = Posets,
      numMutations = 4,
      genotypeMatrix = cbind(1, gMat)
    )
    Results <- ctcbn(bc)
    P <- matrix(0, 219, 24)
    for (i in 1:219) {
        LAMBDA <- as.numeric(Results[[i]]$summary[5:9])
        if (i == 1) {
            DAG <- matrix(0, 0, 0)
        } else {
            DAG <- Posets[[i]]
        }
        P[i, ] <- PathProb_CBN(DAG, LAMBDA, 4)
    }


    ### Step 2: Poset-Level weighting based on the MCMC sampling strategy in BCBN
    Poset_Samples <- bcbn(gMat)

    Poset_Index <- numeric(219)
    for (i in 2:219) {
        d <- dim(Posets[[i]])[1]
        temp_mat <- matrix(0, 4, 4)
        for (j in 1:d) {
            X1 <- Posets[[i]][j, 1]
            X2 <- Posets[[i]][j, 2]
            temp_mat[X1, X2] <- 1
        }
        Poset_Index[i] <- Base2Indexing(temp_mat)
    }

    wB <- numeric(219)
    for (i in 1:100000) {
        Sample_Index <- Base2Indexing(Poset_Samples[[i]])
        INDX <- which(Poset_Index == Sample_Index)
        wB[INDX] <- wB[INDX] + 1
    }

    wB <- wB / sum(wB)


    ### Step 3: Aggregating the Probability Distributions
    PathProb <- apply((wB * P), 2, sum) / sum(wB)

    return(PathProb)
}





#' Jensen_Shannon_Divergence
#'
#' @param Prob1 The first (discrete) probability distribution (vector)
#' @param Prob2 The second (discrete) probability distribution (vector)
#'
#' @return Jensen Shannon Divergence between the two (discrete) probability distributions
#' @export
#'
#' @examples
#' set.seed(100)
#' gMat<-matrix(sample(c(0,1),800,replace = TRUE),200,4)
#' PathCT<-PathProb_Quartet_CTCBN(gMat)
#' PathH<-PathProb_Quartet_HCBN(gMat)
#' JSD<-Jensen_Shannon_Divergence(PathCT,PathH)
Jensen_Shannon_Divergence <- function(Prob1, Prob2) {
    # Prob1: the first probability distribution
    # Prob2: the second probability distribution
    D <- 0
    for (i in 1:length(Prob1)) {
        # if (sum(Prob1[i],na.rm=TRUE)>0){D<-D+Prob1[i]*log2(Prob1[i]/(0.5*Prob1[i]+0.5*Prob2[i]))}
        # if (sum(Prob2[i],na.rm=TRUE)>0){D<-D+Prob2[i]*log2(Prob2[i]/(0.5*Prob1[i]+0.5*Prob2[i]))}
        if (Prob1[i] > 0) {
            D <- D + Prob1[i] * log2(Prob1[i] / (0.5 * Prob1[i] + 0.5 * Prob2[i]))
        }
        if (Prob2[i] > 0) {
            D <- D + Prob2[i] * log2(Prob2[i] / (0.5 * Prob1[i] + 0.5 * Prob2[i]))
        }
    }
    Dv <- (D / 2)
    return(Dv)
}


#' Predictability
#'
#' @param Prob Pathway probability vector
#' @param x The length of genotype vectors
#'
#' @return Predictability
#' @export
#'
#' @examples
#' set.seed(100)
#' gMat<-matrix(sample(c(0,1),800,replace = TRUE),200,4)
#' PathCT<-PathProb_Quartet_CTCBN(gMat)
#' PathH<-PathProb_Quartet_HCBN(gMat)
#' PredC<-Predictability(PathCT,4)
#' PredH<-Predictability(PathH,4)
Predictability <- function(Prob, x) {
    TOT <- 0
    for (i in 1:length(Prob)) {
        if (sum(Prob[i], na.rm = TRUE) > 0) {
            TOT <- TOT - Prob[i] * log(Prob[i])
        }
    }
    Pred <- 1 - (TOT / log(factorial(x))) ## computing the predictability
    return(Pred)
}
