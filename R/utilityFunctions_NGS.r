#' filter input sequences for two patterns
#'
#' @param reads input sequences 
#' @param leftPattern left pattern motive.
#' @param rightPattern right pattern motive.
#' @param maxMismatch maximal number of miss matches.
#' @param prevPatternPos prev pattern posiotion; default is set to NULL.
#' @examples
#' reads <- DNAStringSet(c('ACTGGGTTT','ACCCTGGGTTT'))
#' leftPattern <- 'CT'
#' rightPattern <- 'TTT'
#' maxMismatch <- 0
#' twoPatternReadFilter(reads, leftPattern, rightPattern, maxMismatch)
#' @return list object 
#' @importFrom Biostrings vmatchPattern
#' @export twoPatternReadFilter
twoPatternReadFilter <- function(reads, leftPattern, rightPattern, maxMismatch,
                                 prevPatternPos = NULL){
  vp <- vmatchPattern(leftPattern, reads, max.mismatch = maxMismatch)
  leftStart <- vapply(startIndex(vp), 
                      .dummyFunction, c('startIndex'=0))
  leftEnd <- vapply(endIndex(vp), 
                    .dummyFunction, c('endIndex'=0))
  vp <- vmatchPattern(rightPattern, reads, max.mismatch=maxMismatch)
  rightStart <- vapply(startIndex(vp), 
                       .dummyFunction, c('startIndex'=0))
  rightEnd <- vapply(endIndex(vp), 
                     .dummyFunction, c('endIndex'=0))
  toNA <- which(rightStart < leftEnd)
  rightStart[toNA] <- NA
  patternPositions <- NULL
  if(is.null(prevPatternPos)){
    patternPositions <- cbind(leftStart1 = leftStart, leftEnd1 = leftEnd, 
                              rightStart2 = rightStart, rightEnd2 = rightEnd)
  } else {
    patternPositions <- cbind(leftStart = leftStart, leftEnd = leftEnd, 
                              rightStart = rightStart, rightEnd = rightEnd, 
                              prevPatternPos)
  }
  patternInRead <- !apply(is.na(patternPositions), 1, any)
  reads <- reads[patternInRead]
  patternPositions <- as.data.frame(patternPositions[patternInRead,])
  result <- list(reads = reads, patternPositions = patternPositions)
  return(result)
}

.getRowNameMax <- function(X){
  return(names(sort(X, decreasing = TRUE))[1])
}

.getReadsFromFastq <- function(file){
  f <- FastqFile(file)
  myReads <- sread(readFastq(f))
  close(f)
  return(myReads)
}

.filterByReadLengthMedian <- function(myReads, relMinLength = 0.95, 
                                     relMaxLength = 1.05){
  medianWidth <- median(width(myReads))
  filteredReads <- myReads[which(width(myReads) > medianWidth * relMinLength 
                                 &width(myReads) < medianWidth * relMaxLength)]
  return(filteredReads)
}

.filterReadsWithNs <- function(mySeq){
  readsWithNs <- unique(c(which(grepl('N', vapply(mySeq[['seqFC']], 
                                                  toString, c(FC='')))), 
                          which(grepl('N', vapply(mySeq[['seqNB']], 
                                                  toString, c(NB=''))))))
  if(length(readsWithNs) > 0){
    mySeq[['seqFC']] <- mySeq[['seqFC']][-readsWithNs]
    mySeq[['seqNB']] <- mySeq[['seqNB']][-readsWithNs]
  }
  return(mySeq)
}

.filterForInFrame <- function(mySeq ,minFlycodeLength = 33, 
                             minNanobodyLength = 365){
  width_seqFC <- width(mySeq[['seqFC']])
  width_seqNB <- width(mySeq[['seqNB']])
  mySeq[['seqFC']] <- mySeq[['seqFC']][(width_seqFC %% 3 == 0 
                                        & width_seqFC >= minFlycodeLength 
                                        & width_seqNB %% 3 == 0 
                                        & width_seqNB >= minNanobodyLength)]
  mySeq[['seqNB']] <- mySeq[['seqNB']][(width_seqFC %% 3 == 0 
                                        & width_seqFC >= minFlycodeLength 
                                        & width_seqNB %% 3 == 0 
                                        & width_seqNB >= minNanobodyLength)]
  return(mySeq)
}

.filterForStopCodons <- function(mySeqAS){
  readsWithStops <- unique(c(which(grepl('\\*', 
                                         vapply(mySeqAS[['seqAS_FC']], 
                                                toString, c(FC='')))),
                             which(grepl('\\*', 
                                         vapply(mySeqAS[['seqAS_NB']],
                                                toString, c(NB=''))))))
  mySeqAS[['seqAS_FC']] <- mySeqAS[['seqAS_FC']][-readsWithStops]
  mySeqAS[['seqAS_NB']] <- mySeqAS[['seqAS_NB']][-readsWithStops]
  return(mySeqAS)
}

.filterForStartEndPattern <- function(mySeqAS, myFCPattern = '^GS.*R$'){
  readsWithCorrectFlyCodeEnds <- which(grepl(myFCPattern, 
                                             vapply(mySeqAS[['seqAS_FC']], 
                                                    toString, c(FC=''))))
  mySeqAS[['seqAS_FC']] <- mySeqAS[['seqAS_FC']][readsWithCorrectFlyCodeEnds]
  mySeqAS[['seqAS_NB']] <- mySeqAS[['seqAS_NB']][readsWithCorrectFlyCodeEnds]
  return(mySeqAS)
}

.filterForMinFlycodeFrequency <- function(mySeqAS, minFreq, file, knownFC=c()){
  FCReads <- vapply(mySeqAS[['seqAS_FC']], toString, c(FC=''))
  tmp <- sort(table(FCReads), decreasing = TRUE)
  tmp <- tmp[tmp >= minFreq]
  top_FC <- data.frame(Sequence = names(tmp), Frequency = tmp, 
                       stringsAsFactors = FALSE)
  if(length(knownFC) > 0){
    top_FC[['Name']] <- '-'
    for(j in seq_len(nrow(top_FC))){
      for(k in seq_len(length(knownFC))){
        if(top_FC[['Sequence']][j] == knownFC[[k]])
          top_FC[['Name']][j] <- names(knownFC)[k]
      }
    }
  }
  top_FCFile <- basename(sub('.fastq.gz', '_TopFC.txt', file))
  write.table(top_FC, top_FCFile, sep = '\t', quote = FALSE, row.names = FALSE)
  mySeqAS[['seqAS_FC']] <- AAStringSet(names(tmp))
  names(mySeqAS[['seqAS_FC']]) <- paste('read', 
                                        seq(1,length(mySeqAS[['seqAS_FC']]),1),
                                        sep='_')
  return(mySeqAS)
}

.hit2 <- function(data){
  secondBest <- sort(data, decreasing = TRUE)[2]
  return(secondBest)
}

getConsensusNBs <- function(uniqFCReads, bigTable){
  uniq_FC2NB <- data.frame(Flycode = uniqFCReads, NanoBody = '', 
                           consensusScore = 0, missmatchEvents = 0, 
                          relBestHitFreq = 0, stringsAsFactors = FALSE)
  for (j in seq_len(length(uniqFCReads))){
    NBperFC <- bigTable[bigTable[['Flycode']] == uniqFCReads[j],]
    NBperFC_SS <- AAStringSet(NBperFC[['NanoBody']])
    cMatrix <- consensusMatrix(NBperFC_SS)
    
    hit1 <- apply(cMatrix, 2, max)
    hit2 <- apply(cMatrix, 2, .hit2)
    sum <- apply(cMatrix, 2, sum)
    consensusScore <- mean((hit1 / (hit2 + hit1)))
    missmatches <- sum(sum - hit1)
    identicalNBs <- table(NBperFC[['NanoBody']])
    relBestHitFreq <- max(identicalNBs) / sum(identicalNBs)
    
    if(relBestHitFreq < 0.5){
       mismatchData <- as.numeric(apply(cMatrix, 2, max) < colSums(cMatrix))
       mismatchData <-  gsub('0', '.', as.character(mismatchData))
       mismatchData <- gsub('1', '*', mismatchData)
       mismatchData <- paste(mismatchData, collapse = '')
       fileName <- paste0(names(uniqFCReads)[j], '_', uniqFCReads[j], '.txt')
       write.table(c(NBperFC[['NanoBody']], mismatchData), fileName, sep = '', 
                   quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
    
    uniq_FC2NB[['NanoBody']][j] <- paste(apply(cMatrix, 2, .getRowNameMax), 
                                         collapse = '')
    uniq_FC2NB[['consensusScore']][j] <- consensusScore
    uniq_FC2NB[['missmatchEvents']][j] <- missmatches
    uniq_FC2NB[['relBestHitFreq']][j] <- relBestHitFreq
  }
  return(uniq_FC2NB)
}

.findTopNB <- function(seqAS_NB, n = 100, knownNB){
  top_NB <- head(sort(table(vapply(seqAS_NB, toString, c(NB=''))), 
                      decreasing = TRUE), n = n)
  top_NB <- data.frame(Sequence = names(top_NB), 
                       Frequency = top_NB, stringsAsFactors = FALSE)
  top_NB[['Name']] <- '-'
  for(j in seq_len(nrow(top_NB))){
    for(k in seq_len(length(knownNB))){
      if(top_NB[['Sequence']][j] == knownNB[[k]])
        top_NB[['Name']][j] <- names(knownNB)[k]
    }
  }
  return(top_NB)
}

.dummyFunction <- function(x){if (is.null(x)) NA else x[1]}