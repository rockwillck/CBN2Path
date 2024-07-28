#' CT-CBN
#'
#' @param datasetObj `Spock` object with poset and pattern/lambda data.
#' @param bootstrap_mode Boolean representing bootstrapping mode.
#' @param bootstrap_samples  Number of bootstrap samples (requires `epsilon` > 0, `num_drawn_samples` = 0)
#' @param random_seed Random seed.
#' @param sampling_rate Sampling rate.
#' @param epsilon If between 0 and 1, the fraction of violations allowed per edge. If negative, the interval 0 to 0.5 will be sampled equidistantly with N points.
#' @param num_drawn_samples If > 0, the number of samples to draw from the model. If zero (default), the model will be learned from data.
#' @param num_em_runs Number of em runs.
#' @param outputStem File to output to (defaults to a tempfile).
#'
#' @return A list of output data.
#' @export
#'
#' @examples
#' \dontrun{
#' bc = Spock$new(
#'   poset = read_poset("examples/BC")$sets,
#'   numMutations = read_poset("examples/BC")$mutations,
#'   patternOrLambda = read_pattern("examples/BC")
#' )
#' ctcbn(bc)
#' }
ctcbn <- function(datasetObj,
                  bootstrap_mode = FALSE,
                  bootstrap_samples = 0,
                  random_seed = 0,
                  sampling_rate = 1.0,
                  epsilon = 2,
                  num_drawn_samples = 0,
                  num_em_runs = 1,
                  outputStem = tempfile("output"))
{
  x = .Call(
    "ctcbn_",
    outputStem,
    datasetObj$getPoset(),
    datasetObj$getSecond(num_drawn_samples),
    as.integer(bootstrap_mode),
    as.integer(bootstrap_samples),
    as.integer(random_seed),
    as.double(sampling_rate),
    as.double(epsilon),
    as.integer(num_drawn_samples),
    as.integer(num_em_runs)
  )
  
  splitted = unlist(strsplit(x, "/"))
  
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
  return(outputList)
}