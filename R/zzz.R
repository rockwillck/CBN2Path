#' Get BCBN .tgz path
#'
#' @return internal .tgz path to install rBCBN
#' @export
#'
#' @examples
#' getBCBNinstall()
getBCBNinstall <- function() {
    # Path to the .tgz file within the package
    system.file("extdata", "rBCBN_0.0.0.9000.tgz", package = "CBN2Path")
}
