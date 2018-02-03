#R


#' Compose a peptide with a defined AA sequence frequency
#' @author Christian Panse 
#' @param pool AA distributen.
#' @param cTerm c-Terms
#' @return a AA sequence
#' @export compose_GPGx8cTerm
compose_GPGx8cTerm <- 
  function(pool= c(rep('A', 12), rep('S', 0), rep('T', 12), rep('N', 12), rep('Q', 12), rep('D', 8), 
                   rep('E', 0), rep('V', 12), rep('L', 0), rep('F', 0), rep('Y', 8), rep('W', 0), 
                   rep('G', 12), rep('P', 12)), cTerm=c('VFR','VSR','VFGIR','VSGER')){ 
    paste("GPG", 
          paste(pool[sample(length(pool), 8)], collapse=''), 
          cTerm[sample(length(cTerm), 1)], 
          sep='') 
  }

#' compose a peptide with a defined AA sequence
#' @author Christian Panse 
#' @param pool AA distributen.
#' @param pool AA distributen.
#' @return a AA sequence
#' @export compose_GPx10R
compose_GPx10R <- function(aa_pool1, aa_pool2){ 
  paste("GP", paste(aa_pool1[sample(length(aa_pool1), 2)], collapse=''), 
        paste(aa_pool2[sample(length(aa_pool2), 6)], collapse=''),
        paste(aa_pool1[sample(length(aa_pool1), 2)], collapse=''), "R",
        sep='') 
}

#' plot a LC-MS map
#' @author Christian Panse 
#' @param peptides
#' @return gplots::hist2d 
#' 
.plot_in_silico_LCMS_map <- function(peptides, ...){
  hyd <- unlist(lapply(peptides, function(x){specL::ssrc(x)}))
  pim <- unlist(lapply(peptides, function(x){protViz::parentIonMass(x)}))
  
  pim.range <- range(pim) 
  n.pim <- length(seq(pim.range[1], pim.range[2], by=2))
  n.hyd <- 120
  cm <- c('#44444444', heat.colors(50), 'green')
  n.cm <- length(cm)
  
  h <- gplots::hist2d(hyd, pim,
                      sub = "in-silico lc-ms map (gplots::hist2d)",
                      xlab="hydrophobicity value",
                      ylab='parent ion mass',
                      nbins = c(n.hyd, n.pim),
                      col=cm)
  
  #points(hyd.iRT, pim.iRT, col='cyan')
  legend("topleft", paste("#cells:", n.pim, " x ", n.hyd))
  legend("bottomright", "iRT peptides", col='cyan', pch='o')
  h_counts <- as.numeric(h$counts)
  h_counts <- h_counts[h_counts > 0]
  frequency <- table(h_counts)
  numbers_per_region <- as.numeric(names(frequency))
  
  plot(numbers_per_region, frequency, 
       col = cm[2 + round((n.cm - 1) * (numbers_per_region / max(numbers_per_region)))], 
       xlab = 'number of peptides per region',
       ylab = 'frequency',
       axes=FALSE,
       log='', ...)
  axis(1)
  axis(2, c(1, 10, 100, 1000, 2000), c(1, 10, 100, 1000, 2000))
  abline(v=2.5, col='grey')
  #legend("topleft", paste("#cells:", n.pim, " x ", n.hyd "=", n.pim * n.hyd))
  return(h)
}


get_pool_x8 <- function(){
  aa_pool_x8 <- c(rep('A', 12), rep('S', 0), rep('T', 12), rep('N', 12), rep('Q', 12), rep('D', 8), 
                rep('E', 0), rep('V', 12), rep('L', 0), rep('F', 0), rep('Y', 8), rep('W', 0), 
                rep('G', 12), rep('P', 12))
}


get_pool_1_2_9_10 <- function(){
  aa_pool_1_2_9_10 <- c(rep('A', 8), rep('S', 7), rep('T', 7), rep('N', 6), rep('Q', 6), rep('D', 8), 
                      rep('E', 8), rep('V', 9), rep('L', 6), rep('F', 5), rep('Y', 9), rep('W', 6), 
                      rep('G', 15), rep('P', 0))
}

get_pool_3_8 <- function(){
  aa_pool_3_8 <- c(rep('A', 5), rep('S', 4), rep('T', 5), rep('N', 2), rep('Q', 2), rep('D', 8), 
                 rep('E', 8), rep('V', 7), rep('L', 5), rep('F', 4), rep('Y', 6), rep('W', 4), 
                 rep('G', 12), rep('P', 28))
}


#' Title
#'
#' @return
#' @export
#'
#' @examples
#' @export .getFC
.getFC <- function(){
  
  FC <- read.table(system.file("extdata/FC.tryptic",
                               package = "NestLink"),
                   col.names = c("peptide", "ESP_Prediction"),
                   header = TRUE)
  
  
  
  FC$peptide <- (as.character(FC$peptide))
  idx <- grep ("^GS[ASTNQDEFVLYWGP]{7}(WR|WLTVR|WQEGGR|WLR|WQSR)$", FC$peptide)
  
  FC$cond <- "FC"
  
  FC$pim <- parentIonMass(FC$peptide)
  
  FC <- FC[nchar(FC$peptide) > 2, ]
  
  FC$ssrc <- sapply(FC$peptide, ssrc)
  
  FC$peptideLength <- nchar(as.character(FC$peptide))
  
  unique(FC[idx,])
}


#' Title
#'
#' @return
#' @export
#'
#' @examples
#' export .getNB 
.getNB <- function(){
  
  NB <- read.table(system.file("extdata/NB.tryptic", package = "NestLink"),
                   col.names = c("peptide", "ESP_Prediction"), header = TRUE)
  
  NB$cond <- "NB"
  NB$peptide <- (as.character(NB$peptide))
  NB$pim <- parentIonMass(NB$peptide)
  NB <- NB[nchar(NB$peptide) >2, ]
  NB$ssrc <- sapply(NB$peptide, ssrc)
  NB$peptideLength <- nchar(as.character(NB$peptide))
  # unique(NB)
  NB
}

#' make_it_unambiguous
#'
#' @param x a \code{data.frame} containing a column peptide
#' @return 
#' a \code{data.frame} of unambiguously assignable peptides 
#' (those, which occur only on one nanobody)
#' @export NB.unambiguous
NB.unambiguous <- function(x){
  stopifnot('peptide' %in% names(x))
  
  
  unambiguous <- table(x$peptide) == 1
  
  unambiguous.peptides <- (row.names(unambiguous)[unambiguous])
  
  x$cond <- "NB.unambiguous"
  
  x[x$peptide %in% unambiguous.peptides, ] 
}

#' Title
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
#' @export NB.unique
NB.unique <- function(x){
  stopifnot('peptide' %in% names(x))
  x$cond <- "NB.unique"
  unique(x)
}