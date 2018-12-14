#R

.onAttach <- function(lib, pkg) {
    if (interactive()) {
        version <- packageVersion('NestLink')
        packageStartupMessage("Package 'NestLink' version ", version)
        
        invisible()
    }
}

#' NestLink zzz .onLoad
#'
#' @param libname xx
#' @param pkgname xx
#'
#' @importFrom AnnotationHub query
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom utils read.csv
.onLoad <- function(libname, pkgname) {
    fl <- system.file("extdata", "metadata.csv", package=pkgname)
    metadata <- read.csv(fl, stringsAsFactors=FALSE)

    # eh = ExperimentHub()
    # load(query(eh, c("NestLink", "WU160118.RData"))[[1]])

    # filename <- query(eh, c("NestLink",
    #    "PGexport2_normalizedAgainstSBstandards_Peptides.csv"))[[1]]

    #PGexport2_normalizedAgainstSBstandards_Peptides <- read.csv(filename, header = TRUE, sep=';')
    #query(eh, "NestLink")
    #createHubAccessors(pkgname, titles)
}


#' NestLink getExperimentHubFilename
#'
#' @param filename
#' @return ehubfilename
#' 
#' @export getExperimentHubFilename
#' @importFrom AnnotationHub query
#' @importFrom ExperimentHub ExperimentHub
#' 
#' @examples 
#' (filename <- getExperimentHubFilename("WU160118.RData))
getExperimentHubFilename <- function(filename){
    eh = ExperimentHub()
    suppressMessages({
        filename <- load(query(eh, c("NestLink", "WU160118.RData"))[[1]])

        })
    filename
}