#' CT-CBN Single Batch
#'
#' @param dataset `Spock` object with poset and pattern/lambda data.
#' @param bootstrap_samples  Number of bootstrap samples (requires `epsilon` > 0, `num_drawn_samples` = 0)
#' @param random_seed Random seed.
#' @param sampling_rate Sampling rate.
#' @param epsilon If between 0 and 1, the fraction of violations allowed per edge. If negative, the interval 0 to 0.5 will be sampled equidistantly with N points.
#' @param num_drawn_samples If > 0, the number of samples to draw from the model. If zero (default), the model will be learned from data.
#' @param num_em_runs Number of em runs.
#'
#' @return A list of output data.
#' @export
#'
#' @examples
#' example_path <- getExamples()[1]
#' bc <- Spock$new(
#'     poset = readPoset(example_path)$sets,
#'     numMutations = readPoset(example_path)$mutations,
#'     genotypeMatrix = readPattern(example_path)
#' )
#' ctcbnSingle(bc)
ctcbnSingle <- function(dataset,
                         bootstrap_samples = 0,
                         random_seed = 1,
                         sampling_rate = 1.0,
                         epsilon = 0.0,
                         num_drawn_samples = 0,
                         num_em_runs = 1) {

    if (bootstrap_samples != 0 && !(epsilon > 0 && num_drawn_samples == 0)) {
      stop("bootstrap_samples == 0 requires epsilon > 0 and num_drawn_samples == 0")
    }
    bootstrap_mode = bootstrap_samples != 0

    outputStem <- tempfile("output")
    secondPath <- dataset$getSecond(num_drawn_samples)
    outputs <- list()
    for (i in 1:dataset$getSize()) {
        posetPath <- dataset$getPoset(i)

        x <- .Call(
            "ctcbn_",
            outputStem,
            posetPath,
            secondPath,
            as.integer(bootstrap_mode),
            as.integer(bootstrap_samples),
            as.integer(random_seed),
            as.double(sampling_rate),
            as.double(epsilon),
            as.integer(num_drawn_samples),
            as.integer(num_em_runs)
        )

        splitted <- unlist(strsplit(outputStem, "/"))

        outDir <- paste(splitted[1:(length(splitted) - 1)], collapse = "/")
        outFiles <- (filter_strings_by_start(list.files(outDir), splitted[[length(splitted)]]))

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
            if (endsWith(f, ".time")) {
              outputList$time <- readTime(substring(f, 1, nchar(f) - 5))
            }
            try(file.remove(f), silent = TRUE)
        }

        if (num_drawn_samples == 0) {
          r <- as.numeric(unlist(strsplit(x, " ")))
          labels <- c(c("Poset", "Eps", "Alpha", "Loglike", "lambda_s"), paste0("lambda_", seq(1, length(r) - 5)))
          names(r) <- labels
          outputList$summary <- r
        }

        if (file.exists(paste(posetPath, "poset", sep = "."))) {
            try(file.remove(paste(posetPath, "poset", sep = ".")), silent = TRUE)
        }

        outputs <- append(outputs, list(outputList))
    }
    if (file.exists(paste(secondPath, "lambda", sep = "."))) {
        try(file.remove(paste(secondPath, "lambda", sep = ".")), silent = TRUE)
    }
    if (file.exists(paste(secondPath, "pat", sep = "."))) {
        try(file.remove(paste(secondPath, "pat", sep = ".")), silent = TRUE)
    }

    return(outputs)
}

#' CT-CBN
#'
#' @param datasets Vector of `Spock` objects with poset and pattern/lambda data or a `Spock` object (alias of ctcbnSingle).
#' @param bootstrap_samples  Number of bootstrap samples (requires `epsilon` > 0, `num_drawn_samples` = 0)
#' @param random_seed Random seed.
#' @param sampling_rate Sampling rate.
#' @param epsilon If between 0 and 1, the fraction of violations allowed per edge. If negative, the interval 0 to 0.5 will be sampled equidistantly with N points.
#' @param num_drawn_samples If > 0, the number of samples to draw from the model. If zero (default), the model will be learned from data.
#' @param num_em_runs Number of em runs.
#' @param n_cores Maximum number of threads to use to parallelize.
#'
#' @return A matrix of results.
#' @export
#'
#' @examples
#' example_path <- getExamples()[1]
#' bc <- Spock$new(
#'     poset = readPoset(example_path)$sets,
#'     numMutations = readPoset(example_path)$mutations,
#'     genotypeMatrix = readPattern(example_path)
#' )
#' ctcbn(bc)
ctcbn <- function(datasets,
                  bootstrap_samples = 0,
                  random_seed = 1,
                  sampling_rate = 1.0,
                  epsilon = 2,
                  num_drawn_samples = 0,
                  num_em_runs = 1,
                  n_cores = 1) {
    if (length(datasets) < n_cores) {
        message(paste("Number of datasets was less than number of cores. Using number of datasets (", length(datasets), ") as thread count.", sep = ""))
    }
    registerDoMC(cores = min(length(datasets), n_cores))

    if (inherits(datasets, "Spock")) {
        return(ctcbnSingle(datasets, bootstrap_samples, random_seed, sampling_rate, epsilon, num_drawn_samples, num_em_runs))
    }
    output_stems <- replicate(length(datasets), tempfile("output"))
    rowLength <- -1
    done <- 0
    outMatrixBuf <- vector("list", length(datasets))
    dataI <- 1
    rets <- foreach(dataI = 1:length(datasets)) %dopar% {
        out <- ctcbnSingle(datasets[[dataI]], bootstrap_samples, random_seed, sampling_rate, epsilon, num_drawn_samples, num_em_runs)
        list(i = dataI, row = out)
    }

    return(rets)
}
