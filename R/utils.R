#' Read a .pat file
#'
#' @param filestem The filename of the .pat file without the .pat suffix.
#'
#' @return A matrix.
#' @export
#'
#' @examples
#' bcPath = get_examples()[1]
#' read_pattern(bcPath)
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
#' bcPath = get_examples()[1]
#' read_poset(bcPath)
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
#' bcPath = get_examples()[1]
#' read_lambda(bcPath)
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

#' Filters strings by start
#' @param strings A list of strings.
#' @param start_substring Starting string to filter by.
#' @noRd
filter_strings_by_start <- function(strings, start_substring) {
  filtered_strings <- grep(paste0("^", start_substring), strings, value = TRUE)
  return(filtered_strings)
}

#' Get paths to examples
#'
#' @return A vector of paths
#' @export
#'
#' @examples
#' get_examples()
get_examples <- function() {
  examples = c("BC", "CRC", "hiv", "HSD", "prostate", "RCC")
  gsub(".poset", "", lapply(examples, function(x) system.file("extdata", paste(x, ".poset", sep=""), package="rCBN")))
}

pad_list <- function(list, length) {
  vec = unlist(list)
  length_out = length - length(vec)
  if (length_out > 0) {
    c(vec, rep(NA, length_out))
  } else {
    vec
  }
}
