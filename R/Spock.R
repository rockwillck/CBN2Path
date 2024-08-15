#' Poset and pattern/lambda data
#'
#' @description
#' A data class containing poset and pattern/lambda matrices.
#'
#' @details
#' Use the read_ methods to feed data from files.
#' @export
Spock = R6::R6Class("Spock", list(
  
  #' @field poset Poset matrix.
  poset = matrix(),
  #' @field numMutations Number of mutations.
  numMutations = 0,
  #' @field patternOrLambda Pattern/lambda matrix.
  patternOrLambda = matrix(),
  
  #' @description
  #' Create a new Spock object.
  #' @param poset Poset matrix.
  #' @param numMutations Number of mutations.
  #' @param patternOrLambda Pattern matrix or lambda matrix.
  #' @return A new `Spock` object.
  initialize = function (poset, numMutations, patternOrLambda) {
    stopifnot(is.matrix(poset))
    stopifnot(is.numeric(numMutations))
    stopifnot(is.matrix(patternOrLambda))
    self$poset = poset
    self$numMutations = numMutations
    self$patternOrLambda = patternOrLambda
  },
  
  #' @description
  #' Write poset data to a tempfile.
  #' @return File path to tempfile.
  getPoset = function() {
    output = ""
    if (nrow(self$poset) == 2) {
      for (i in 1:ncol(self$poset)) {
        output = paste(output, paste(self$poset[i, 1], self$poset[i, 2]), sep = "")
        output = paste(output, "\n")
      }
    }
    fileC = (paste(paste(
      as.character(self$numMutations), output, sep = "\n"
    ), "0", sep = ""))
    return(temp_file(ext = ".poset", fileC = fileC))
  },
  
  #' @description
  #' Write pattern/lambda data to a tempfile.
  #' @param n Number of drawn samples.
  #' @return File path to tempfile.
  getSecond = function(n) {
    if (n == 0) {
      return(self$getPattern())
    } else {
      return(self$getLambda())
    }
  },
  
  #' @description
  #' Write pattern data to a tempfile.
  #' @return File path to tempfile.
  getPattern = function() {
    fileC = (paste(
      paste(as.character(nrow(self$patternOrLambda)), as.character(ncol(self$patternOrLambda)), sep = " "),
      matrix_to_string(self$patternOrLambda),
      sep = "\n"
    ))
    
    return(temp_file(ext = ".pat", fileC = fileC))
  },
  
  #' @description
  #' Write lambda data to a tempfile.
  #' @return File path to tempfile.
  getLambda = function() {
    output = ""
    for (i in 1:length(self$patternOrLambda)) {
      output = paste(output, self$patternOrLambda[[i]], sep = "")
      output = paste(output, "\n")
    }
    fileC = (output)
    return(temp_file(ext = ".lambda", fileC = fileC))
  }
))