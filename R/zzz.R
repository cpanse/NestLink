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


#' getExperimentHubFilename
#'
#' @param filename of the aws s3 blob.
#' @return the file name of the local ExperimentHub.
#' 
#' @export getExperimentHubFilename
#' @importFrom AnnotationHub query
#' @importFrom ExperimentHub ExperimentHub
#' 
#' @examples 
#' fl <- system.file("extdata", "metadata.csv", package="NestLink")   
#' metadata <- read.csv(fl, stringsAsFactors=FALSE)
#'      metadata$Title     
#'      
#' lapply(metadata$RDataPath, getExperimentHubFilename)
getExperimentHubFilename <- function(filename){
    suppressMessages({
        eh = ExperimentHub()
    
        filename <- query(eh, c("NestLink", filename))[[1]]
        })
    filename
}