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
#' example_path = get_examples()[1]
#' bc = Spock$new(
#'   poset = read_poset(example_path)$sets,
#'   numMutations = read_poset(example_path)$mutations,
#'   genotypeMatrix = read_pattern(example_path)
#' )
#' hcbn(bc)
hcbn <- function(datasetObj,
                 anneal = FALSE,
                 temp = 0,
                 annealing_steps = 0,
                 outputStem = tempfile("output"))
{
  posetPath = datasetObj$getPoset()
  secondPath = datasetObj$getSecond(0)
  x = .Call(
    "hcbn_",
    outputStem,
    posetPath,
    secondPath,
    as.integer(anneal),
    as.double(temp),
    as.integer(annealing_steps)
  )
  
  splitted = unlist(strsplit(x, "/"))
  
  outDir = paste(splitted[1:(length(splitted) - 1)], collapse = "/")
  outFiles = (filter_strings_by_start(list.files(outDir), splitted[[length(splitted)]]))
  
  outputList = list()
  for (f in outFiles) {
    f = paste(c(outDir, f), collapse = "/")
    if (endsWith(f, ".poset")) {
      outputList$poset = read_poset(substring(f, 1, nchar(f) - 6))
    }
    if (endsWith(f, ".pat")) {
      outputList$pattern = read_pattern(substring(f, 1, nchar(f) - 4))
    }
    if (endsWith(f, ".lambda")) {
      outputList$lambda = read_lambda(substring(f, 1, nchar(f) - 7))
    }
    file.remove(f)
  }
  
  if (file.exists(paste(posetPath,"poset",sep="."))) {
    file.remove(paste(posetPath,"poset",sep="."))
  }
  if (file.exists(paste(secondPath,"lambda",sep="."))) {
    file.remove(paste(secondPath,"lambda",sep="."))
  }
  if (file.exists(paste(secondPath,"pat",sep="."))) {
    file.remove(paste(secondPath,"pat",sep="."))
  }
  
  return(outputList)
}
