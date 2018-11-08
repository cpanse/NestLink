#R

#' write FASTA
#'
#' @param x a \code{nanobodyFlycodeLinking} class computed by \code{\link{runNGSAnalysis}}.

#' @param ... just passed
#'
#' @return sprintf stream
#' @export nanobodyFlycodeLinking.as.fasta
#'
#' @author Lennart Opitz, Christian Panse 2018
#' @examples
#' \dontrun{
#' FF <- c('~/__projects/2018/20181030-p1875/NL7-1.extendedFrags_uniqNB2FC.txt', 
#' '~/__projects/2018/20181030-p1875/NL7-2.extendedFrags_uniqNB2FC.txt')
#' 
#' X <- do.call('rbind', lapply(FF, read.table, header=T, sep='\t'))
#' cat(uniqNB2FC.as.fasta(X), file=tempfile(fileext=".fasta"))
#' }
nanobodyFlycodeLinking.as.fasta <- function(x, ...){
  idx <- seq(1, nrow(x))
  sprintf(">NB%04d_NL7Idx1 FC%d %s\n%s\n", idx, x$FlycodeCount[idx],
          x$NB[idx], gsub(",", "", x$AssociatedFlycodes[idx]))
}


#' Object Summaries of S3 class \code{nanobodyFlycodeLinking}
#'
#' @param object a \code{nanobodyFlycodeLinking} class computed by \code{\link{runNGSAnalysis}}.
#'
#' @return a data.frame object
#' @export nanobodyFlycodeLinking.summary
nanobodyFlycodeLinking.summary <- function(object){
  
  # number of NB
  # number FC
  # number of AAs
  # cat(length(object))
  # cat()
  cat ("to be implemented.")
}

