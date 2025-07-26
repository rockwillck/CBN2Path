#' Poset and pattern/lambda data
#'
#' @description A data class containing poset and pattern/lambda matrices.
#'
#' @details Use the read_ methods to feed data from files.
#' @return a Spock object
#' @export
#'
#' @examples
#' example_path <- get_examples()[1]
#' bc <- Spock$new(
#'     poset = read_poset(example_path)$sets,
#'     numMutations = read_poset(example_path)$mutations,
#'     genotypeMatrix = read_pattern(example_path)
#' )
Spock <- R6::R6Class("Spock", list(

    #' @field poset Poset matrix.
    poset = NULL,
    #' @field numMutations Number of mutations.
    numMutations = 0,
    #' @field genotypeMatrix Genotype matrix.
    genotypeMatrix = matrix(),
    #' @field lambda Lambda list.
    lambda = NULL,

    #' @description Create a new Spock object.
    #' @param poset Poset matrix or list of poset matrices.
    #' @param numMutations Number of mutations.
    #' @param genotypeMatrix Genotype matrix.
    #' @param lambda Lambda list.
    #' @return A new `Spock` object.
    initialize = function(poset, numMutations, genotypeMatrix, lambda = NULL) {
        stopifnot(is.matrix(poset) || is.list(poset))
        stopifnot(is.numeric(numMutations))
        stopifnot(is.matrix(genotypeMatrix))
        if (is.matrix(poset)) {
            self$poset <- list(poset)
        } else {
            self$poset <- poset
        }
        self$numMutations <- numMutations
        self$genotypeMatrix <- genotypeMatrix
        self$lambda <- lambda
    },

    #' @description Get the number of posets.
    #' @return Number of posets.
    getSize = function() {
        return(length(self$poset))
    },

    #' @description Write poset data to a tempfile.
    #' @param index Index of poset.
    #' @return File path to tempfile.
    getPoset = function(index = 1) {
        output <- ""
        pos <- self$poset[[index]]
        if (ncol(pos) == 2) {
            for (i in 1:nrow(pos)) {
                output <- paste(output, paste(pos[i, 1], pos[i, 2]), sep = "")
                output <- paste(output, "\n")
            }
        }
        fileC <- (paste(paste(
            as.character(self$numMutations), output,
            sep = "\n"
        ), "0", sep = ""))
        return(temp_file(ext = ".poset", fileC = fileC))
    },

    #' @description Write pattern/lambda data to a tempfile.
    #' @param n Number of drawn samples.
    #' @return File path to tempfile.
    getSecond = function(n) {
        if (n == 0) {
            return(self$getPattern())
        } else {
            return(self$getLambda())
        }
    },

    #' @description Write pattern data to a tempfile.
    #' @return File path to tempfile.
    getPattern = function() {
        fileC <- (paste(
            paste(as.character(nrow(self$genotypeMatrix)), as.character(ncol(self$genotypeMatrix)), sep = " "),
            matrix_to_string(self$genotypeMatrix),
            sep = "\n"
        ))

        return(temp_file(ext = ".pat", fileC = fileC))
    },

    #' @description Write lambda data to a tempfile.
    #' @return File path to tempfile.
    getLambda = function() {
        if (is.null(self$lambda)) {
            output <- ""
            for (i in 1:length(self$genotypeMatrix)) {
                output <- paste(output, self$genotypeMatrix[[i]], sep = "")
                output <- paste(output, "\n")
            }
            fileC <- (output)
        } else {
            output <- ""
            for (i in 1:length(self$lambda)) {
                output <- paste(output, self$lambda[[i]], sep = "")
                output <- paste(output, "\n")
            }
            fileC <- (output)
        }
        return(temp_file(ext = ".lambda", fileC = fileC))
    }
))
