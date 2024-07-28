#' Read a .pat file
#'
#' @param filestem The filename of the .pat file without the .pat suffix.
#'
#' @return A matrix.
#' @export
#'
#' @examples
#' \dontrun{
#' read_pattern("examples/BC")
#' }
read_pattern <- function(filestem) {
  lines <- readLines(paste(suppressWarnings(normalizePath(filestem)),
                           ".pat",
                           sep = ""))

  # Parse dimensions from the first line
  dimensions <- as.numeric(strsplit(lines[1], "\\s+")[[1]][2:3])
  num_rows <- dimensions[1]
  num_cols <- dimensions[2]

  # Read matrix data from subsequent lines
  matrix_data <- sapply(lines[2:length(lines)], function(line)
    as.numeric(strsplit(line, "\\s+")[[1]]))

  # Convert matrix data into a matrix
  matrix_data <- matrix(matrix_data, nrow = num_rows, byrow = FALSE)
  return(t(matrix_data))
}

#' Read a .poset file
#'
#' @param filestem The filename of the .poset file without the .poset suffix.
#'
#' @return A list containing the number of mutations and a matrix.
#' @export
#'
#' @examples
#' \dontrun{
#' read_poset("examples/BC")
#' }
read_poset <- function(filestem) {
  lines <- readLines(paste(suppressWarnings(normalizePath(filestem)), ".poset", sep =
                             ""))
  if (trimws(lines[[2]]) == "0") {
    return(list(mutations = as.numeric(lines[[1]]), sets = matrix()))
  }
  allSets = c()
  for (i in 3:length(lines) - 1) {
    allSets = c(allSets, as.numeric(unlist(
      strsplit(lines[[i]], " ")
    )[[1]]))
    allSets = c(allSets, as.numeric(unlist(
      strsplit(lines[[i]], " ")
    )[[2]]))
  }

  return(list(mutations = as.numeric(lines[[1]]), sets = matrix(allSets, ncol=2, byrow = TRUE)))
}

#' Read a .lambda file
#'
#' @param filestem The filename of the .lambda file without the .lambda suffix.
#'
#' @return A matrix.
#' @export
#'
#' @examples
#' \dontrun{
#' read_lambda("examples/BC")
#' }
read_lambda <- function(filestem) {
  lines <- unlist(as.numeric(readLines(
    paste(suppressWarnings(normalizePath(filestem)), ".lambda", sep = "")
  )))
  return(matrix(lines, ncol = 1))
}

#' Writes a matrix to a string
#' @param mat A matrix.
#' @noRd
matrix_to_string <- function(mat) {
  # Get the dimensions of the matrix
  nrows <- nrow(mat)
  ncols <- ncol(mat)

  # Initialize an empty string to store the result
  result <- ""

  # Loop through each row
  for (i in 1:nrows) {
    res <- paste(mat[i, ], collapse = " ")
    # Loop through each column
    # for (j in 1:ncols) {
    #   # Append the matrix element to the result string
    #   result <- paste(result, mat[i, j], sep = " ")
    # }
    # Add a newline character at the end of each row (except the last row)
    if (i == 1) {
      result <- res

    } else if (i <= nrows) {
      result <- paste(result, res, sep = "\n")
    }
  }

  return(result)
}

#' Write to a tempfile
#' @param ext File extension.
#' @param fileC File contents.
#' @noRd
temp_file = function(ext, fileC) {
  tf = tempfile("tempo", fileext = ext)[[1]]
  tfObj = file(tf)
  write(fileC, tfObj)
  close(tfObj)
  return(substr(tf, 0, nchar(tf) - nchar(ext)))
}

#' Poset and pattern/lambda data
#'
#' @description
#' A data class containing poset and pattern/lambda matrices.
#'
#' @details
#' Use the read_ methods to feed data from files.
#' @export
Spock = R6::R6Class("Spock", list(

  #' @field poset Poset list.
  poset = list(),
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

#' Filters strings by start
#' @param strings A list of strings.
#' @param start_substring Starting string to filter by.
#' @noRd
filter_strings_by_start <- function(strings, start_substring) {
  filtered_strings <- grep(paste0("^", start_substring), strings, value = TRUE)
  return(filtered_strings)
}
