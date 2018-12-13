#R

.onAttach <- function(lib, pkg) {
    if (interactive()) {
        version <- packageVersion('NestLink')
        packageStartupMessage("Package 'NestLink' version ", version)
        
        invisible()
    }
}

#' NestLinkzzz
#'
#' @param libname xx
#' @param pkgname xx
#'
#' @return createHubAccessors
#' @import ExperimentHub
#' @importFrom utils read.csv
.onLoad <- function(libname, pkgname) {
    fl <- system.file("extdata", "metadata.csv", package=pkgname)
    metadata <- read.csv(fl, stringsAsFactors=FALSE)
    #library(ExperimentHub)
    #eh = ExperimentHub()
    #query(eh, "NestLink")
    #createHubAccessors(pkgname, titles)
}
