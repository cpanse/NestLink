#R


FF <- c('~/__projects/2018/20181030-p1875/NL7-1.extendedFrags_uniqNB2FC.txt', 
        '~/__projects/2018/20181030-p1875/NL7-2.extendedFrags_uniqNB2FC.txt')



#' write FASTA
#'
#' @param x 
#' @param file 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' FF <- c('~/__projects/2018/20181030-p1875/NL7-1.extendedFrags_uniqNB2FC.txt', 
#' '~/__projects/2018/20181030-p1875/NL7-2.extendedFrags_uniqNB2FC.txt')
#' 
#' X<-do.call('rbind', lapply(FF, read.table, header=T, sep='\t'))
#' cat(xxx.as.fasta(X), file="/tmp/XX.fasta")
#' }
xxx.as.fasta <- function(x, file = tempfile(fileext=".fasta"), ...){
  idx<-1:nrow(x)
  sprintf(">NB%04d_NL7Idx1 FC%d %s\n%s\n", idx, X$FlycodeCount[idx], X$NB[idx], gsub(",", "", X$AssociatedFlycodes[idx]))
}
