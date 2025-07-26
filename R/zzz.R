#' Get BCBN .tgz path
#'
#' @return internal .tgz path to install rBCBN
#' @export
#'
#' @examples
#' getBCBNinstall()
getBCBNinstall <- function() {
    files = list.files(system.file("extdata", package = "CBN2Path"), full.names = TRUE)
    match = files[grepl("^rBCBN", basename(files))]
    
    match[1]
}

#' Get needed BCBN version
#' @NoRd
getBCBNVersion <- function() {
  pathsplit = strsplit(getBCBNinstall(),"/")[[1]]
  packageVersionR = pathsplit[[length(pathsplit)]]
  pkgVn = substr(packageVersionR, 7, nchar(packageVersionR)-4)
  
  pkgVn
}

#' Check if installed BCBN is correct
#' @NoRd
checkBCBNVersion <- function() {
  identical(as.character(utils::packageVersion("rBCBN")), getBCBNVersion())
}