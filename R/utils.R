#' Read a .pat file
#'
#' @param fileStem The filename of the .pat file without the .pat suffix.
#'
#' @return A matrix.
#' @export
#'
#' @examples
#' bcPath <- getExamples()[1]
#' readPattern(bcPath)
readPattern <- function(fileStem) {
    lines <- suppressWarnings(readLines(paste(suppressWarnings(normalizePath(fileStem)),
        ".pat",
        sep = ""
    )))

    # Parse dimensions from the first line
    dimensions <- as.numeric(strsplit(lines[1], "\\s+")[[1]][2:3])
    numRows <- dimensions[1]
    numCols <- dimensions[2]

    # Read matrix data from subsequent lines
    matrixData <- sapply(lines[2:length(lines)], function(line) {
        as.numeric(strsplit(line, "\\s+")[[1]])
    })

    # Convert matrix data into a matrix
    matrixData <- matrix(matrixData, nrow = numRows, byrow = FALSE)
    return(t(matrixData))
}

#' Read a .time file
#'
#' @param fileStem The filename of the .time file without the .time suffix.
#'
#' @return A matrix.
#' @export
#'
#' @examples
#' bcPath <- getExamples()[1]
#' readPattern(bcPath)
readTime <- function(fileStem) {
  lines <- suppressWarnings(readLines(paste(suppressWarnings(normalizePath(fileStem)),
                           ".time",
                           sep = ""
  )))

  # Parse dimensions from the first line
  dimensions <- as.numeric(strsplit(lines[1], "\\s+")[[1]][2:3])
  numRows <- dimensions[1]
  numCols <- dimensions[2]

  # Read matrix data from subsequent lines
  matrixData <- sapply(lines[2:length(lines)], function(line) {
    as.numeric(strsplit(line, "\\s+")[[1]])
  })

  # Convert matrix data into a matrix
  matrixData <- matrix(matrixData, nrow = numRows, byrow = FALSE)
  return(t(matrixData))
}

#' Read a .poset file
#'
#' @param fileStem The filename of the .poset file without the .poset suffix.
#'
#' @return A list containing the number of mutations and a matrix.
#' @export
#'
#' @examples
#' bcPath <- getExamples()[1]
#' readPoset(bcPath)
readPoset <- function(fileStem) {
    lines <- suppressWarnings(readLines(paste(suppressWarnings(normalizePath(fileStem)),
                                              ".poset",
                                              sep = ""
    )))
    if (trimws(lines[[2]]) == "0") {
        return(list(mutations = as.numeric(lines[[1]]), sets = matrix()))
    }
    allSets <- c()
    for (i in 3:length(lines) - 1) {
        allSets <- c(allSets, as.numeric(unlist(
            strsplit(lines[[i]], " ")
        )[[1]]))
        allSets <- c(allSets, as.numeric(unlist(
            strsplit(lines[[i]], " ")
        )[[2]]))
    }

    return(list(mutations = as.numeric(lines[[1]]), sets = matrix(allSets, ncol = 2, byrow = TRUE)))
}

#' Read a .lambda file
#'
#' @param fileStem The filename of the .lambda file without the .lambda suffix.
#'
#' @return A matrix.
#' @export
#'
#' @examples
#' bcPath <- getExamples()[1]
#' readLambda(bcPath)
readLambda <- function(fileStem) {
    lines <- suppressWarnings(
      unlist(as.numeric(readLines(
        paste(normalizePath(fileStem), ".lambda", sep = "")
    )))
    )
    return(matrix(lines, ncol = 1))
}

#' Writes a matrix to a string
#' @param mat A matrix.
#' @noRd
matrixToString <- function(mat) {
    # Get the dimensions of the matrix
    nRows <- nrow(mat)
    nCols <- ncol(mat)

    # Initialize an empty string to store the result
    result <- ""

    # Loop through each row
    for (i in 1:nRows) {
        res <- paste(mat[i, ], collapse = " ")
        # Add a newline character at the end of each row (except the last row)
        if (i == 1) {
            result <- res
        } else if (i <= nRows) {
            result <- paste(result, res, sep = "\n")
        }
    }

    return(result)
}

#' Write to a tempfile
#' @param ext File extension.
#' @param fileC File contents.
#' @noRd
tempFile <- function(ext, fileC) {
    tf <- tempfile("tempo", fileext = ext)[[1]]
    tfObj <- file(tf)
    write(fileC, tfObj)
    close(tfObj)
    return(substr(tf, 0, nchar(tf) - nchar(ext)))
}

#' Filters strings by start
#' @param strings A list of strings.
#' @param startSubstring Starting string to filter by.
#' @noRd
filterStringsByStart <- function(strings, startSubstring) {
    filteredStrings <- grep(paste0("^", startSubstring), strings, value = TRUE)
    return(filteredStrings)
}

#' Get paths to examples
#'
#' @return A vector of paths
#' @export
#'
#' @examples
#' getExamples()
getExamples <- function() {
    examples <- c("BC", "CRC", "hiv", "HSD", "prostate", "RCC")
    gsub(".poset", "", lapply(examples, function(x) system.file("extdata", paste(x, ".poset", sep = ""), package = "CBN2Path")))
}

padList <- function(list, length) {
    vec <- unlist(list)
    lengthOut <- length - length(vec)
    if (lengthOut > 0) {
        c(vec, rep(NA, lengthOut))
    } else {
        vec
    }
}

getParents <- function(poset, i){
  allParents <- which(poset[, i] == 1)
  sort(allParents)
}

matrixPower <- function (X, n){
  ## unfortunately R doesn't have a matrix power function like MATLAB.
  ## This function is taken from packages Biodem and is originally called mtx.exp.
  if (n != round(n)) {
    n <- round(n)
    warning("rounding exponent `n' to", n)
  }
  phi <- diag(nrow = nrow(X))
  pot <- X
  while (n > 0) {
    if (n%%2)
      phi <- phi %*% pot
    n <- n%/%2
    pot <- pot %*% pot
  }
  return(phi)
}

