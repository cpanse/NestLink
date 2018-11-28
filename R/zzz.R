#R

.onAttach <- function(lib, pkg) {
    if (interactive()) {
        version <- packageVersion('NestLink')
        packageStartupMessage("Package 'NestLink' version ", version)
        
        invisible()
    }
}

.onLoad <- function(libname, pkgname) {
    fl <- system.file("extdata", "metadata.csv", package=pkgname)
    titles <- read.csv(fl, stringsAsFactors=FALSE)$Title
    createHubAccessors(pkgname, titles)
}
