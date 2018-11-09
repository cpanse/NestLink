#R


#' Compose a FlyCode GSx7cTerm Amino Acid Sequence
#'
#' @description composes, out of a given input distributen,
#' a random sampled amino acid sequence.
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
  function(pool=c(rep('A', 18), rep('S', 6), rep('T', 12),
    rep('N', 1),
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
    rep('Q', 12), rep('D', 8), rep('E', 0), rep('V', 12), rep('L', 0),
    rep('F', 0), rep('Y', 8), rep('W', 0), rep('G', 12), rep('P', 12)),
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
#' aa_pool_1_2_9_10 <- c(rep('A', 8), rep('S', 7), rep('T', 7), rep('N', 6), 
#' rep('Q', 6), rep('D', 8), rep('E', 8), rep('V', 9), rep('L', 6), rep('F', 5),
#' rep('Y', 9), rep('W', 6), rep('G', 15), rep('P', 0))
#' 
#' aa_pool_3_8 <- c(rep('A', 5), rep('S', 4), rep('T', 5), rep('N', 2),
#' rep('Q', 2), rep('D', 8), rep('E', 8), rep('V', 7), rep('L', 5), rep('F', 4),
#' rep('Y', 6), rep('W', 4), rep('G', 12), rep('P', 28))
#' 
#' compose_GPx10R(aa_pool_1_2_9_10, aa_pool_3_8)
#' (FlyCodes <- replicate(10, compose_GPx10R(aa_pool_1_2_9_10, aa_pool_3_8)))
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


#' plot a LC-MS map of a given set of amino acid sequences
#' @author Christian Panse 
#' @param peptides a vector of pepitdes.
#' @param ... pass through the plot method.
#' @importFrom graphics abline axis barplot legend plot
#' @importFrom grDevices dev.off heat.colors png
#' @importFrom gplots hist2d
#' @details TODO(cp): consider using hexbin using ggplot2
#' ggplot facet_wrap aes geom_point 
#' @importFrom protViz parentIonMass ssrc
#' @export plot_in_silico_LCMS_map
#' @examples 
#' set.seed(1)
#' par(mfrow=c(2,1));
#' FlyCodes <- replicate(10000, compose_GPGx8cTerm())
#' rv <- plot_in_silico_LCMS_map(FlyCodes)
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
       col = cm[2 + round((n.cm - 1) * 
         (numbers_per_region / max(numbers_per_region)))], 
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


#' Read FlyCodes (FC)
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
getFC <- function(
  pattern = "^GS[ASTNQDEFVLYWGP]{7}(WR|WLTVR|WQEGGR|WLR|WQSR)$"){
  
  FC <- read.table(system.file("extdata/FC.tryptic",
                               package = "NestLink"),
                   col.names = c("peptide", "ESP_Prediction"),
                   header = TRUE)
  
  
  
  FC$peptide <- (as.character(FC$peptide))
  idx <- grep (pattern, FC$peptide)
  
  FC$cond <- "FC"
  
  FC$pim <- parentIonMass(FC$peptide)
  
  FC <- FC[nchar(FC$peptide) > 2, ]
  
  FC$ssrc <- ssrc(FC$peptide)
  # FC$ssrc <- vapply(FC$peptide, FUN=ssrc, FUN.VALUE = ssrc("ELVISR"))
  FC$peptideLength <- nchar(as.character(FC$peptide))
  
  unique(FC[idx,])
}



#' Read NanoBodies (NB)
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
  NB$ssrc <- ssrc(NB$peptide)
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
#' class(WU160118)
#' PATTERN <- "^GS[ASTNQDEFVLYWGP]{7}(WR|WLTVR|WQEGGR|WLR|WQSR)$"
#' idx <- grepl(PATTERN, WU160118$pep_seq)
#' WU <- WU160118[idx & WU160118$pep_score > 25,]
#' 
#' library(lattice)
#' histogram(~RTINSECONDS| datfilename, data = WU, type='count')
NULL


#' F255744 Mascot Search results
#' 
#' @name F255744
#' @docType data
#' @author Pascal Egloff \email{p.egloff@imm.uzh.ch}
#' @seealso \href{https://fgcz-bfabric.uzh.ch/bfabric/userlab/show-resource.html?id=409912}{F255744}
#' @keywords data
#' @examples 
#' class(F255744)
#' hist(F255744$RTINSECONDS)
#' hist(F255744$RTINSECONDS[F255744$pep_score > 20])
NULL


#' PGexport results
#' 
#' @name PGexport
#' @docType data
#' @source 
#' \url{https://fgcz-bfabric.uzh.ch}
#' \itemize{ 
#' \item{Workunit : 158716 - QEXACTIVEHF_1
#'  20170919_16_62465_nl5idx1-3_6titratecoli.raw
#'  20170919_05_62465_nl5idx1-3_6titratecoli.raw}
#' \item{Workunit : 158717 - QEXACTIVEHF_1
#'  20170919_14_62466_nl5idx1-3_7titratesmeg.raw
#'  20170919_09_62466_nl5idx1-3_7titratesmeg.raw}}
#' @author Pascal Egloff \email{p.egloff@imm.uzh.ch}
#' @keywords data
#' @examples 
#' filename <- system.file(
#'   "extdata/PGexport2_normalizedAgainstSBstandards_Peptides.csv",
#'   package = "NestLink")
#' P <- read.csv(filename, header = TRUE, sep=';')
#' P <- P[P$Modifications == '', ]
#' P <- P[,c('Accession', 'Sequence', 
#' "X20170919_05_62465_nl5idx1.3_6titratecoli", 
#' "X20170919_16_62465_nl5idx1.3_6titratecoli",  
#' "X20170919_09_62466_nl5idx1.3_7titratesmeg", 
#' "X20170919_14_62466_nl5idx1.3_7titratesmeg")]
#' names(P)<-c('Accession','Sequence','coli1', 'coli2', 'smeg1', 'smeg2')
#' P<- P[grep("^P[0-9][A-Z][0-9]", P$Accession), ] 
#' 
#' P$FCset_ng <- NA
#' P$FCset_ng[P$Accession %in% c('P1A4', 'P1B4', 'P1C4',
#'   'P1D4', 'P1E4', 'P1F4')] <- 92
#' P$FCset_ng[P$Accession %in% c('P1A5', 'P1B5', 'P1C5',
#'   'P1D5', 'P1G4', 'P1H4')] <- 295
#' P$FCset_ng[P$Accession %in% c('P1A6', 'P1B6', 'P1E5',
#'   'P1F5', 'P1G5', 'P1H5')] <- 943
#' P$FCset_ng[P$Accession %in% c('P1C6', 'P1D6', 'P1E6', 
#'   'P1F6', 'P1G6', 'P1H6')] <- 3017
#'   
#'  P$coli1 <- (log(P$coli1,2) - mean(log(P$coli1,2))) / sd(log(P$coli1,2))
#'  P$coli2 <- (log(P$coli2,2) - mean(log(P$coli2,2))) / sd(log(P$coli2,2))
#'  P$smeg1 <- (log(P$smeg1,2) - mean(log(P$smeg1,2))) / sd(log(P$smeg1,2))
#'  P$smeg2 <- (log(P$smeg2,2) - mean(log(P$smeg2,2))) / sd(log(P$smeg2,2))
#'  
#'   O<-P
#'   b <- boxplot(df<-cbind(P$coli1 - P$coli2, P$coli1 - P$smeg1, 
#'   P$coli1 - P$smeg2,P$coli2 - P$smeg1, P$coli2 - P$smeg2 , P$smeg1 - P$smeg2),
#'    ylab='normalized log2ratios', ylim = c(-1,1), axes=FALSE,
#'     main=paste("ConcGr = all"))
#' axis(1, 1:6, c('coli[12]', 'coli1-smeg1', 'coli1-smeg2', 'coli2-smeg1',
#' 'coli2- smeg2','smeg[12]'))
#' abline(h=0, col='red')
#' box()
#' axis(2)
#' axis(3, 1:6, b$n)
#' outliers.idx <- sapply(1:length(b$group), function(i){
#'   q <- df[, b$group[i]] == b$out[i];
#'   text(b$group[i], b$out[i], P[q, 2], pos=4, cex=0.4);
#'    text(b$group[i], b$out[i], P[q, 1], pos=2, cex=0.4);
#'    which(q)}
#'    )
NULL