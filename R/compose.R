#R


#' Compose a FlyCode GSx7cTerm Amino Acid Sequence
#'
#' @description composes, out of a given input distributen, a random sampled amino acid sequence.
#' @param pool a vector of amino acids.
#' @param cTerm a vector of a sequence suffix.
#'
#' @return a amino acid sequence, e.g., GSAPTTVFGWLTVR.
#' @export compose_GSx7cTerm
#' 
#' @examples  
#'
#'  sample.size <- 100
#'  #
#'  ## Compose a GSXXXXXXX(WR|WLTVR|WQGGER|WQSR|WLR) peptide
#'  set.seed(2)
#'  FC.GSx7cTerm <- replicate(sample.size, compose_GSx7cTerm())
#'  ## Some Sanity Checks
#'  table(FC.GSx7cTerm)
#'  stopifnot(length(FC.GSx7cTerm) == 100)
#'  FC.PATTERN <- "^GS[ASTNQDEFVLYWGP]{7}(WR|WLTVR|WQEGGR|WLR|WQSR)$"
#'  stopifnot(
#'    length(FC.GSx7cTerm[grepl(FC.PATTERN, FC.GSx7cTerm)]) 
#'      == sample.size)
#' 
#' @author Christian Panse <cp@fgcz.ethz.ch> 2015
#' @export compose_GSx7cTerm
compose_GSx7cTerm <- 
  function(pool=c(rep('A', 18), rep('S', 6), rep('T', 12), rep('N', 1),
                  rep('Q', 1), rep('D', 11), rep('E', 11), rep('V', 12), 
                  rep('L', 2), rep('F', 1), rep('Y', 4), rep('W', 1), 
                  rep('G', 8), rep('P', 12)), 
           cTerm=c('WR','WLTVR','WQEGGR','WQSR','WLR')){ 
    
    paste("GS", 
          paste(pool[sample(length(pool), 7)], collapse=''), 
          cTerm[sample(length(cTerm), 1)], 
          sep='') 
  }


#' Compose a peptide with a defined AA sequence frequency
#' @author Christian Panse, 2015
#' @param pool AA distributen.
#' @param cTerm c-Terms
#' @author Christian Panse <cp@fgcz.ethz.ch> 2015
#' @return a AA sequence
#' @examples 
#' set.seed(1)
#' compose_GPGx8cTerm()
#' (FlyCodes <- replicate(10, compose_GPGx8cTerm()))
#' plot(parentIonMass(FlyCodes) ~ssrc(FlyCodes))
#' @export compose_GPGx8cTerm
compose_GPGx8cTerm <- 
  function(pool= c(rep('A', 12), rep('S', 0), rep('T', 12), rep('N', 12),
                   rep('Q', 12), rep('D', 8), rep('E', 0), rep('V', 12),
                   rep('L', 0), rep('F', 0), rep('Y', 8), rep('W', 0), 
                   rep('G', 12), rep('P', 12)), 
           cTerm=c('VFR','VSR','VFGIR','VSGER')){ 
    paste("GPG", 
          paste(pool[sample(length(pool), 8)], collapse=''), 
          cTerm[sample(length(cTerm), 1)], 
          sep='') 
  }

#' Compose a peptide with a defined AA sequence
#' @author Christian Panse <cp@fgcz.ethz.ch> 2015 
#' @param aa_pool1 AA distributen.
#' @param aa_pool2 AA distributen.
#' @examples 
#' set.seed(1)
#' compose_GPx10R()
#' (FlyCodes <- replicate(10, compose_GPx10R()))
#' plot(parentIonMass(FlyCodes) ~ssrc(FlyCodes))
#' @return a AA sequence
#' @export compose_GPx10R
compose_GPx10R <- function(aa_pool1, aa_pool2){ 
  paste("GP", paste(aa_pool1[sample(length(aa_pool1), 2)], collapse=''), 
        paste(aa_pool2[sample(length(aa_pool2), 6)], collapse=''),
        paste(aa_pool1[sample(length(aa_pool1), 2)], collapse=''), "R",
        sep='') 
}

.plot_in_silico_LCMS_map <- function(peptides, ...){
  plot_in_silico_LCMS_map(peptides, ...)
}


#' plot a LC-MS map
#' @author Christian Panse 
#' @param peptides a vector of pepitdes.
#' @param ... pass through the plot method.
#' @importFrom graphics abline axis barplot legend plot
#' @importFrom grDevices dev.off heat.colors png
#' @importFrom gplots hist2d
# ggplot2 ggplot facet_wrap aes geom_point 
#' @importFrom protViz parentIonMass ssrc
#' @examples 
#' set.seed(1)
#' plot_in_silico_LCMS_map(FlyCodes <- replicate(10, compose_GPx10R()))
#' @return gplots::hist2d a gplot 2d histogram
plot_in_silico_LCMS_map <- function(peptides, ...){
  hyd <- unlist(lapply(peptides, function(x){ssrc(x)}))
  pim <- unlist(lapply(peptides, function(x){parentIonMass(x)}))
  
  pim.range <- range(pim) 
  n.pim <- length(seq(pim.range[1], pim.range[2], by=2))
  n.hyd <- 120
  cm <- c('#44444444', heat.colors(50), 'green')
  n.cm <- length(cm)
  
  h <- hist2d(hyd, pim,
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
  aa_pool_x8 <- c(rep('A', 12), rep('S', 0), rep('T', 12), rep('N', 12),
                  rep('Q', 12), rep('D', 8), rep('E', 0), rep('V', 12),
                  rep('L', 0), rep('F', 0), rep('Y', 8), rep('W', 0), 
                rep('G', 12), rep('P', 12))
}


get_pool_1_2_9_10 <- function(){
  aa_pool_1_2_9_10 <- c(rep('A', 8), rep('S', 7), rep('T', 7), rep('N', 6),
                        rep('Q', 6), rep('D', 8),  rep('E', 8), rep('V', 9),
                        rep('L', 6), rep('F', 5), rep('Y', 9), rep('W', 6), 
                      rep('G', 15), rep('P', 0))
}

get_pool_3_8 <- function(){
  aa_pool_3_8 <- c(rep('A', 5), rep('S', 4), rep('T', 5), rep('N', 2),
                   rep('Q', 2), rep('D', 8), rep('E', 8), rep('V', 7), 
                   rep('L', 5), rep('F', 4), rep('Y', 6), rep('W', 4), 
                 rep('G', 12), rep('P', 28))
}


#' Read FlyCodes
#'
#' @return a data.frame of FlyCodes
#' @param pattern a regular expression FlyCode pattern 
#' @author Christian Panse <cp@fgcz.ethz.ch> 2015
#' @importFrom protViz parentIonMass
#' @importFrom utils read.table write.table
#' @examples 
#' FC <- getFC()
#' dim(FC)
#' 
#' @export getFC
getFC <- function(pattern = "^GS[ASTNQDEFVLYWGP]{7}(WR|WLTVR|WQEGGR|WLR|WQSR)$"){
  
  FC <- read.table(system.file("extdata/FC.tryptic",
                               package = "NestLink"),
                   col.names = c("peptide", "ESP_Prediction"),
                   header = TRUE)
  
  
  
  FC$peptide <- (as.character(FC$peptide))
  idx <- grep (pattern, FC$peptide)
  
  FC$cond <- "FC"
  
  FC$pim <- parentIonMass(FC$peptide)
  
  FC <- FC[nchar(FC$peptide) > 2, ]
  
  FC$ssrc <- sapply(FC$peptide, ssrc)
  
  FC$peptideLength <- nchar(as.character(FC$peptide))
  
  unique(FC[idx,])
}



#' Read NanoBodies
#'
#' @return a data.frame of NBs
#' @examples
#' NB <- getNB()
#' dim(NB)
#' @export getNB 
getNB <- function(){
  
  NB <- read.table(system.file("extdata/NB.tryptic", package = "NestLink"),
                   col.names = c("peptide", "ESP_Prediction"), header = TRUE)
  
  NB$cond <- "NB"
  NB$peptide <- (as.character(NB$peptide))
  NB$pim <- parentIonMass(NB$peptide)
  NB <- NB[nchar(NB$peptide) >2, ]
  NB$ssrc <- sapply(NB$peptide, ssrc)
  NB$peptideLength <- nchar(as.character(NB$peptide))
  NB
}

#' Determine unambiguous NBs
#'
#' @param x a \code{data.frame} containing a column peptide
#' @return a data.frame 
#' a \code{data.frame} of unambiguously assignable peptides 
#' (those, which occur only on one nanobody)
#' @examples 
#' NB <- getNB()
#' dim(NB.unambiguous(NB))
#' @export NB.unambiguous
#' @importFrom grDevices pdf
#' @importFrom stats aggregate median
#' @importFrom utils packageVersion
NB.unambiguous <- function(x){
  stopifnot('peptide' %in% names(x))
  
  
  unambiguous <- table(x$peptide) == 1
  
  unambiguous.peptides <- (row.names(unambiguous)[unambiguous])
  
  x$cond <- "NB.unambiguous"
  
  x[x$peptide %in% unambiguous.peptides, ] 
}

#' make NB table unique
#'
#' @param x a data.frame
#' @importFrom protViz ssrc
#' @return a data.frame
#' @export NB.unique
#'
#' @examples
#' NB <- getNB()
#' dim(NB.unique(NB))
#' @export NB.unique
NB.unique <- function(x){
  stopifnot('peptide' %in% names(x))
  x$cond <- "NB.unique"
  unique(x)
}

#' WU160118 Mascot Search results
#' 
#' @name WU160118
#' @docType data
#' @author Christian Panse 
#' @references \url{https://fgcz-bfabric.uzh.ch/bfabric/userlab/show-workunit.html?id=160118}
#' @keywords data
#' @examples 
#' data(WU160118)
#' PATTERN <- "^GS[ASTNQDEFVLYWGP]{7}(WR|WLTVR|WQEGGR|WLR|WQSR)$"
#' idx <- grepl(PATTERN, WU160118$pep_seq)
#' WU <- WU160118[idx & WU160118$pep_score > 25,]
NULL


#' F255744 Mascot Search results
#' 
#' @name F255744
#' @docType data
#' @author Pascal Egloff \email{p.egloff@imm.uzh.ch}
#' @references \url{http://fgcz-mascot-server.uzh.ch/mascot/cgi/master_results_2.pl?file=..%2Fdata%2F20170819%2FF255744.dat}
#' @keywords data
NULL
