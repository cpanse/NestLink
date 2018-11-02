#R



#' write FASTA
#'
#' @param x 
#' @param file 
#' @param ... 
#'
#' @return
#' @export nanobodyFlycodeLinking.as.fasta
#'
#' @author Lennart Opitz, Christian Panse 2018
#' @examples
#' \dontrun{
#' FF <- c('~/__projects/2018/20181030-p1875/NL7-1.extendedFrags_uniqNB2FC.txt', 
#' '~/__projects/2018/20181030-p1875/NL7-2.extendedFrags_uniqNB2FC.txt')
#' 
#' X<-do.call('rbind', lapply(FF, read.table, header=T, sep='\t'))
#' cat(uniqNB2FC.as.fasta(X), file="/tmp/XX.fasta")
#' }
nanobodyFlycodeLinking.as.fasta <- function(x, file = tempfile(fileext=".fasta"), ...){
  idx <- 1:nrow(x)
  message(paste("writting to file ", file))
  sprintf(">NB%04d_NL7Idx1 FC%d %s\n%s\n", idx, x$FlycodeCount[idx], x$NB[idx], gsub(",", "", x$AssociatedFlycodes[idx]))
}
