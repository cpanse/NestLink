#R

.onAttach <- function(lib, pkg) {
    if (interactive()) {
        version <- packageVersion('NestLink')
        packageStartupMessage("Package 'NestLink' version ", version)
        
        invisible()
    }
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
