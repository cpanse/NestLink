#R

#' Write FASTA
#'
#' @param x a \code{nanobodyFlycodeLinking} S3 object computed by
#' \code{\link{runNGSAnalysis}}.
#' @param file a filename
#' @param ... just passed
#' @return sprintf stream
#' @export nanobodyFlycodeLinking.as.fasta
#'
#' @author Lennart Opitz, Christian Panse 2018
#' @examples
#' library(ExperimentHub)
#' eh <- ExperimentHub()
#' f <- query(eh, c("NestLink", "nanobodyFlycodeLinkage.RData"))[[1]]
#' load(f)
#' summary(nanobodyFlycodeLinkage.sample)
#' nanobodyFlycodeLinking.as.fasta(nanobodyFlycodeLinkage.sample)
nanobodyFlycodeLinking.as.fasta <- function(x, file = NULL, ...) {
    if (!is.nanobodyFlycodeLinking(x)) {
        warning("object is not of class nanobodyFlycodeLinking")
    }
    idx <- seq(1, nrow(x))
    
    fasta <- sprintf(
        ">NB%04d FC%d %s\n%s\n",
        idx,
        x$FlycodeCount[idx],
        x$NB[idx],
        gsub(",", "", x$AssociatedFlycodes[idx])
    )
    
    if (!is.null(file)) {
        cat(fasta, file, sep = '')
        message(paste("FASTA written to", file))
    } else{
        return(fasta)
    }
}

is.nanobodyFlycodeLinking <- function(object) {
    sum(object$FlycodeCount) == sum(vapply(strsplit(object$AssociatedFlycodes,
                                                    ","), length, 1))
}

#' Object Summaries of S3 class \code{nanobodyFlycodeLinking}
#'
#' @param object a \code{nanobodyFlycodeLinking} class computed
#' by \code{\link{runNGSAnalysis}}.
#'
#' @return a data.frame object
#' @export nanobodyFlycodeLinking.summary
#' @examples
#' library(ExperimentHub)
#' eh <- ExperimentHub()
#' f <- query(eh, c("NestLink", "nanobodyFlycodeLinkage.RData"))[[1]]
#' load(f)
#' summary(nanobodyFlycodeLinkage.sample)
nanobodyFlycodeLinking.summary <- function(object) {
    if (!is.nanobodyFlycodeLinking(object)) {
        warning("object is not of class nanobodyFlycodeLinking")
    }
    df <- data.frame(
        "number.Nanobodies" = nrow(object),
        "number.Flycodes" = sum(object$FlycodeCount),
        "number.AminoAcids" = sum(nchar(
            gsub(",", "", object$AssociatedFlycodes)
        ))
    )
    df
}
