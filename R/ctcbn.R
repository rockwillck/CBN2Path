#' CT-CBN Single Batch
#'
#' @param dataset `Spock` object with poset and pattern/lambda data.
#' @param bootstrap_mode Boolean representing bootstrapping mode.
#' @param bootstrap_samples  Number of bootstrap samples (requires `epsilon` > 0, `num_drawn_samples` = 0)
#' @param random_seed Random seed.
#' @param sampling_rate Sampling rate.
#' @param epsilon If between 0 and 1, the fraction of violations allowed per edge. If negative, the interval 0 to 0.5 will be sampled equidistantly with N points.
#' @param num_drawn_samples If > 0, the number of samples to draw from the model. If zero (default), the model will be learned from data.
#' @param num_em_runs Number of em runs.
#' @param output_stem File to output to (defaults to a tempfile).
#'
#' @return A list of output data.
#' @export
#'
#' @examples
#' example_path = get_examples()[1]
#' bc = Spock$new(
#'   poset = read_poset(example_path)$sets,
#'   numMutations = read_poset(example_path)$mutations,
#'   patternOrLambda = read_pattern(example_path)
#' )
#' ctcbn_single(bc)
ctcbn_single <- function(dataset,
                  bootstrap_mode = FALSE,
                  bootstrap_samples = 0,
                  random_seed = 0,
                  sampling_rate = 1.0,
                  epsilon = 2,
                  num_drawn_samples = 0,
                  num_em_runs = 1,
                  output_stem = tempfile("output"))
{
  x = .Call(
    "ctcbn_",
    output_stem,
    dataset$getPoset(),
    dataset$getSecond(num_drawn_samples),
    as.integer(bootstrap_mode),
    as.integer(bootstrap_samples),
    as.integer(random_seed),
    as.double(sampling_rate),
    as.double(epsilon),
    as.integer(num_drawn_samples),
    as.integer(num_em_runs)
  )
  
  splitted = unlist(strsplit(output_stem, "/"))
  
  outDir = paste(splitted[1:(length(splitted) - 1)], collapse = "/")
  outFiles = (filter_strings_by_start(list.files(outDir), splitted[[length(splitted)]]))
  
  outputList = list()
  for (file in outFiles) {
    file = paste(c(outDir, file), collapse = "/")
    if (endsWith(file, ".poset")) {
      outputList$poset = read_poset(substring(file, 1, nchar(file) - 6))
    }
    if (endsWith(file, ".pat")) {
      outputList$pattern = read_pattern(substring(file, 1, nchar(file) - 4))
    }
    if (endsWith(file, ".lambda")) {
      outputList$lambda = read_lambda(substring(file, 1, nchar(file) - 7))
    }
  }
  r = as.numeric(unlist(strsplit(x, " ")))
  labels = c(c("Poset", "Eps", "Alpha", "Loglike", "lambda_s"), paste0("lambda_", seq(1, length(r) - 5)))
  names(r) = labels
  outputList$row = r
  return(outputList)
}

#' CT-CBN
#'
#' @param datasets Vector of `Spock` objects with poset and pattern/lambda data or a `Spock` object (alias of ctcbn_single).
#' @param bootstrap_mode Boolean representing bootstrapping mode.
#' @param bootstrap_samples  Number of bootstrap samples (requires `epsilon` > 0, `num_drawn_samples` = 0)
#' @param random_seed Random seed.
#' @param sampling_rate Sampling rate.
#' @param epsilon If between 0 and 1, the fraction of violations allowed per edge. If negative, the interval 0 to 0.5 will be sampled equidistantly with N points.
#' @param num_drawn_samples If > 0, the number of samples to draw from the model. If zero (default), the model will be learned from data.
#' @param num_em_runs Number of em runs.
#' @param output_stems Files to output to (defaults to tempfiles).
#' @param n_cores Maximum number of threads to use to parallelize.
#'
#' @return A matrix of results.
#' @export
#'
#' @examples
#' example_path = get_examples()[1]
#' bc = Spock$new(
#'   poset = read_poset(example_path)$sets,
#'   numMutations = read_poset(example_path)$mutations,
#'   patternOrLambda = read_pattern(example_path)
#' )
#' ctcbn(bc)
#' ctcbn(c(bc, bc, bc))
ctcbn <- function(datasets,
                    bootstrap_mode = FALSE,
                    bootstrap_samples = 0,
                    random_seed = 0,
                    sampling_rate = 1.0,
                    epsilon = 2,
                    num_drawn_samples = 0,
                    num_em_runs = 1,
                    output_stems = NULL,
                    n_cores = 1) {
  
  if (length(datasets) < n_cores) {
    message(paste("Number of datasets was less than number of cores. Using number of datasets (", length(datasets), ") as thread count.", sep = ""))
  }
  registerDoMC(cores=min(length(datasets), n_cores))
  
  if (inherits(datasets, "Spock")) {
    return(ctcbn_single(datasets, bootstrap_mode, bootstrap_samples, random_seed, sampling_rate, epsilon, num_drawn_samples, num_em_runs, ifelse(is.null(output_stems), tempfile("output"), output_stems[[1]])))
  }
  if (is.null(output_stems)) {
    output_stems = replicate(length(datasets), tempfile("output"))
  }
  rowLength = -1
  done = 0
  outMatrixBuf = vector("list", length(datasets))
  dataI = 1
  rets = foreach(dataI = 1:length(datasets)) %dopar% {
    out = ctcbn_single(datasets[[dataI]], bootstrap_mode, bootstrap_samples, random_seed, sampling_rate, epsilon, num_drawn_samples, num_em_runs, output_stems[[dataI]])
    r = unname(out$row)
    list(i = dataI, row = list(r))
  }
  
  rowLengths = sapply(rets, function(x) length(unlist(x$row)))
  rowLength = max(rowLengths)
  if (length(unique(rowLengths)) > 1) {
    warning("Result vectors are not the same length. Padding with NA.")
  }
  outMatrix = matrix(nrow = length(datasets), ncol = rowLength)
  for (row in rets) {
    outMatrix[row$i,] = pad_list(row$row, rowLength)
  }
  outMatrix
}
