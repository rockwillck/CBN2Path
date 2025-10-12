#' H-CBN Single Batch
#'
#' @param datasetObj `Spock` object with poset and pattern/lambda data.
#' @param anneal If `TRUE`, performes a simulated annealing run starting from the poset
#' @param temp Temperature of simulated annealing.
#' @param annealingSteps Number of simulated annealing steps.
#' @param epsilon Value of eps for CT-CBN model selection. Requires both pattern and lambda data in input `Spock`.
#'
#' @return A list of output data.
#' @export
#'
#' @examples
#' examplePath <- getExamples()[1]
#' bc <- Spock$new(
#'     poset = readPoset(examplePath)$sets,
#'     numMutations = readPoset(examplePath)$mutations,
#'     genotypeMatrix = readPattern(examplePath)
#' )
#' hcbnSingle(bc)
hcbnSingle <- function(datasetObj,
                        anneal = FALSE,
                        temp = 0,
                        annealingSteps = 0,
                        epsilon = 2) {
    if ((epsilon > 0 && epsilon != 2) && is.null(datasetObj$lambda)) {
        stop("Spock object should have lambda list if non-zero epsilon provided.")
    }
    outputStem <- tempfile("output")
    secondPath <- datasetObj$getSecond(0)
    thirdPath <- datasetObj$getSecond(1)
    outputs <- list()
    for (i in 1:datasetObj$getSize()) {
        posetPath <- datasetObj$getPoset(i)

        x <- .Call(
            "hcbn_",
            outputStem,
            posetPath,
            secondPath,
            thirdPath,
            as.integer(anneal),
            as.double(temp),
            as.integer(annealingSteps),
            as.double(epsilon)
        )

        splitted <- unlist(strsplit(outputStem, "/"))

        outDir <- paste(splitted[1:(length(splitted) - 1)], collapse = "/")
        outFiles <- (filterStringsByStart(list.files(outDir), splitted[[length(splitted)]]))

        outputList <- list()
        for (f in outFiles) {
            f <- paste(c(outDir, f), collapse = "/")
            if (endsWith(f, ".poset")) {
                outputList$poset <- readPoset(substring(f, 1, nchar(f) - 6))
            }
            if (endsWith(f, ".pat")) {
                outputList$pattern <- readPattern(substring(f, 1, nchar(f) - 4))
            }
            if (endsWith(f, ".lambda")) {
                outputList$lambda <- readLambda(substring(f, 1, nchar(f) - 7))
            }
            # suppressWarnings(file.remove(f))
        }

        if (file.exists(paste(posetPath, "poset", sep = "."))) {
            # suppressWarnings(file.remove(paste(posetPath, "poset", sep = ".")))
        }

        r <- as.numeric(unlist(strsplit(x, " ")))
        labels <- c(c("Poset", "Eps", "Alpha", "Loglike", "lambda_s"), paste0("lambda_", seq(1, length(r) - 5)))
        names(r) <- labels
        outputList$summary <- r

        outputs <- append(outputs, list(outputList))
    }
    if (file.exists(paste(thirdPath, "lambda", sep = "."))) {
        # suppressWarnings(file.remove(paste(thirdPath, "lambda", sep = ".")))
    }
    if (file.exists(paste(secondPath, "pat", sep = "."))) {
        # suppressWarnings(file.remove(paste(secondPath, "pat", sep = ".")))
    }

    if (length(outputs) == 1) {
      return(outputs[[1]])
    }

    return(outputs)
}

#' H-CBN
#'
#' @param datasets Vector of `Spock` objects with poset and pattern/lambda data or a `Spock` object (alias of hcbnSingle).
#' @param anneal If `TRUE`, performes a simulated annealing run starting from the poset
#' @param temp Temperature of simulated annealing.
#' @param annealingSteps Number of simulated annealing steps.
#' @param epsilon Value of eps for CT-CBN model selection. Requires both pattern and lambda data in input `Spock`.
#' @param nCores Maximum number of threads to use to parallelize.
#'
#' @return A matrix of results.
#' @export
#'
#' @examples
#' examplePath <- getExamples()[1]
#' bc <- Spock$new(
#'     poset = readPoset(examplePath)$sets,
#'     numMutations = readPoset(examplePath)$mutations,
#'     genotypeMatrix = readPattern(examplePath)
#' )
#' hcbn(bc)
#' hcbn(c(bc, bc, bc))
hcbn <- function(datasets,
                 anneal = FALSE,
                 temp = 0,
                 annealingSteps = 0,
                 epsilon = 2,
                 nCores = 1) {

    if (inherits(datasets, "Spock") && length(datasets$poset) == 1) {
      return(hcbnSingle(datasets, anneal, temp, annealingSteps, epsilon))
    } else if (inherits(datasets, "Spock")) {
      datasets = sapply(datasets$poset, \(poset) Spock$new(
        poset = poset,
        numMutations = datasets$numMutations,
        genotypeMatrix = datasets$genotypeMatrix
      ))
    }

  if (length(datasets) < nCores) {
    message(paste("Number of datasets was less than number of cores. Using number of datasets (", length(datasets), ") as thread count.", sep = ""))
  }

    outputStems <- replicate(length(datasets), tempfile("output"))
    rowLength <- -1
    done <- 0
    outMatrixBuf <- vector("list", length(datasets))

    if (exists("MulticoreParam", mode = "function")) {
      p <- MulticoreParam(workers = min(length(datasets), nCores))
      rets <- bplapply(datasets, \(x) hcbnSingle(x, anneal, temp, annealingSteps, epsilon), BPOPTIONS = bpoptions(progressbar = TRUE), BPPARAM = p)
    } else {
      message("MulticoreParam not found — running sequentially.")
      rets <- lapply(datasets, \(x) hcbnSingle(x, anneal, temp, annealingSteps, epsilon))
    }
    
    return(rets)
}
