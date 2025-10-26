#' CT-CBN Single Batch
#'
#' @param dataset `Spock` object with poset and pattern/lambda data.
#' @param bootstrapSamples Number of bootstrap samples (requires `epsilon` > 0, `numDrawnSamples` = 0)
#' @param randomSeed Random seed.
#' @param samplingRate Sampling rate.
#' @param epsilon If between 0 and 1, the fraction of violations allowed per edge. If negative, the interval 0 to 0.5 will be sampled equidistantly with N points and the output `Spock` will include multiple resulting posets.
#' @param numDrawnSamples If > 0, the number of samples to draw from the model. If zero (default), the model will be learned from data.
#' @param numEmRuns Number of em runs.
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
#' ctcbnSingle(bc)
ctcbnSingle <- function(dataset,
                         bootstrapSamples = 0,
                         randomSeed = 1,
                         samplingRate = 1.0,
                         epsilon = 2,
                         numDrawnSamples = 0,
                         numEmRuns = 1) {

    if (bootstrapSamples != 0 && !(epsilon > 0 && numDrawnSamples == 0)) {
      stop("bootstrapSamples == 0 requires epsilon > 0 and numDrawnSamples == 0")
    }
    bootstrapMode = bootstrapSamples != 0

    outputStem <- tempfile("output")
    secondPath <- dataset$getSecond(numDrawnSamples)
    outputs <- list()
    for (i in 1:dataset$getSize()) {
        posetPath <- dataset$getPoset(i)

        x <- .Call(
            "ctcbn_",
            outputStem,
            posetPath,
            secondPath,
            as.integer(bootstrapMode),
            as.integer(bootstrapSamples),
            as.integer(randomSeed),
            as.double(samplingRate),
            as.double(epsilon),
            as.integer(numDrawnSamples),
            as.integer(numEmRuns)
        )

        splitted <- unlist(strsplit(outputStem, "/"))

        outDir <- paste(splitted[1:(length(splitted) - 1)], collapse = "/")
        outFiles <- (filterStringsByStart(list.files(outDir), splitted[[length(splitted)]]))

        outputList <- list()
        outputList$poset <- list()
        # print(outFiles)
        for (f in outFiles) {
            f <- paste(c(outDir, f), collapse = "/")
            if (endsWith(f, ".poset")) {
                outputList$poset <- c(outputList$poset, list(readPoset(substring(f, 1, nchar(f) - 6))))
            }
            if (endsWith(f, ".pat")) {
                outputList$pattern <- readPattern(substring(f, 1, nchar(f) - 4))
            }
            if (endsWith(f, ".lambda")) {
                outputList$lambda <- readLambda(substring(f, 1, nchar(f) - 7))
            }
            if (endsWith(f, ".time")) {
              outputList$time <- readTime(substring(f, 1, nchar(f) - 5))
            }

            if (endsWith(f, ".summary")) {
              x <- readLines(f)
            }
            # suppressWarnings(file.remove(f))
        }

        if (length(outputList$poset) == 1) {
          outputList$poset = outputList$poset[[1]]
        }

        if (numDrawnSamples == 0) {

          r <- lapply(x, \(row) as.numeric(unlist(strsplit(row, " "))))
          labels <- c(c("Poset", "Eps", "Alpha", "Loglike", "lambda_s"), paste0("lambda_", seq(1, length(r[[1]]) - 5)))

          df <- as.data.frame(do.call(rbind, r))
          colnames(df) <- labels

          outputList$summary <- df
        }

        if (file.exists(paste(posetPath, "poset", sep = "."))) {
            # suppressWarnings(file.remove(paste(posetPath, "poset", sep = ".")))
        }

        outputs <- append(outputs, list(outputList))
    }
    if (file.exists(paste(secondPath, "lambda", sep = "."))) {
        # suppressWarnings(file.remove(paste(secondPath, "lambda", sep = ".")))
    }
    if (file.exists(paste(secondPath, "pat", sep = "."))) {
        # suppressWarnings(file.remove(paste(secondPath, "pat", sep = ".")))
    }

    if (length(outputs) == 1) {
      return(outputs[[1]])
    }

    return(outputs)
}

#' CT-CBN
#'
#' @param datasets Vector of `Spock` objects with poset and pattern/lambda data or a `Spock` object (alias of ctcbnSingle).
#' @param bootstrapSamples Number of bootstrap samples (requires `epsilon` > 0, `numDrawnSamples` = 0)
#' @param randomSeed Random seed.
#' @param samplingRate Sampling rate.
#' @param epsilon If between 0 and 1, the fraction of violations allowed per edge. If negative, the interval 0 to 0.5 will be sampled equidistantly with N points and the output `Spock` will include multiple resulting posets.
#' @param numDrawnSamples If > 0, the number of samples to draw from the model. If zero (default), the model will be learned from data.
#' @param numEmRuns Number of em runs.
#' @param nCores Maximum number of threads to use to parallelize.
#' @param progressBar Print out progress bar; default is FALSE
#'
#' @return A matrix of results.
#' @export
#'
#' @examples
#' examplePath <- getExamples()[3]
#' bc <- Spock$new(
#'     poset = readPoset(examplePath)$sets,
#'     numMutations = readPoset(examplePath)$mutations,
#'     genotypeMatrix = readPattern(examplePath)
#' )
#' ctcbn(bc)
#' ctcbn(c(bc, bc, bc))
ctcbn <- function(datasets,
                  bootstrapSamples = 0,
                  randomSeed = 1,
                  samplingRate = 1.0,
                  epsilon = 2,
                  numDrawnSamples = 0,
                  numEmRuns = 1,
                  nCores = 1,
                  progressBar = FALSE) {
    if (inherits(datasets, "Spock") && length(datasets$poset) == 1) {
      return(ctcbnSingle(datasets, bootstrapSamples, randomSeed, samplingRate, epsilon, numDrawnSamples, numEmRuns))
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


    if (exists("MulticoreParam", mode = "function") || exists("SnowParam", mode = "function")) {
      if(Sys.info()["sysname"] == "Windows") {
        p <- SnowParam(workers = min(length(datasets), nCores))
      } else {
        p <- MulticoreParam(workers = min(length(datasets), nCores))
      }
      rets <- bplapply(datasets, \(x) ctcbnSingle(x, bootstrapSamples, randomSeed, samplingRate, epsilon, numDrawnSamples, numEmRuns), BPOPTIONS = bpoptions(progressbar = progressBar), BPPARAM = p)
    } else {
      message("Parallelization not found - running sequentially.")
      rets <- lapply(datasets, \(x) ctcbnSingle(x, bootstrapSamples, randomSeed, samplingRate, epsilon, numDrawnSamples, numEmRuns))
    }

    return(rets)
}
