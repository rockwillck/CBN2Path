#' H-CBN
#'
#' @param datasetObj `Spock` object with poset and pattern/lambda data.
#' @param anneal If `TRUE`, performes a simulated annealing run starting from the poset
#' @param temp Temperature of simulated annealing.
#' @param annealing_steps Number of simulated annealing steps.
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
#' hcbn(bc)
#' }
hcbn <- function(datasetObj,
                 anneal = FALSE,
                 temp = 0,
                 annealing_steps = 0,
                 outputStem = tempfile("output"))
{
  x = .Call(
    "hcbn_",
    outputStem,
    datasetObj$getPoset(),
    datasetObj$getSecond(0),
    as.integer(anneal),
    as.double(temp),
    as.integer(annealing_steps)
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