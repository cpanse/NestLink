#' NGS linkage workflow
#' @description  performs the NGS filtering workflow to get high quality
#' FlyCode and Nanobody sequences linkage.
#' @param file sequence file path
#' @param param list of input parameters, explained in details paragraph below.
#' @author Lennart Opitz <lopitz@fgcz.ethz.ch>, 2019
#' @details 
#' The elements of the parameter list object is described as follows:
#' \itemize{
#' \item{\code{NB_Linker1}} {nucleotide sequence of the linker left to the nanobody.}
#' \item{\code{NB_Linker2}} {nucleotide sequence of the linker right to the nanobody.}
#' \item{\code{ProteaseSite}} {nucleotide sequence left to the flycode.}
#' \item{\code{FC_Linker}} {nucleotide sequence right to the flycode.}
#' \item{\code{knownNB}} {known nanobody sequences in the experiment.}
#' \item{\code{nReads}} {number of Reads from the start of fastq file to process.}
#' \item{\code{minRelBestHitFreq}} {minimal fraction of the dominant nanobody for a specific flycode.}
#' \item{\code{minConsensusScore}} {minimal fraction per sequence position in nanabody consensus sequence calculation.}
#' \item{\code{maxMismatch}} {number of accepted mismatches for all pattern search steps.}
#' \item{\code{minNanobodyLength}} {minimal nanobody length in [nt].}
#' \item{\code{minFlycodeLength}} { minimal flycode length in [nt].}
#' \item{\code{FCminFreq}} {minimal number of subreads for a specific flycode to keep it in the analysis.}
#' }
#' missing elements are replace by the example provided values.
#' 
#' @return uniqNB2FC dataframe
#' @export runNGSAnalysis
#' @import ShortRead
#' @import Biostrings
#' @examples
#' library(ExperimentHub)
#' eh <- ExperimentHub()
#' expFile <- query(eh, c("NestLink", "NL42_100K.fastq.gz"))[[1]]
#' knownNB_File <- query(eh, c("NestLink", "knownNB.txt"))[[1]]    
#' knownNB_data <- read.table(knownNB_File, sep='\t', header = TRUE,
#'     row.names = 1, stringsAsFactors = FALSE)
#' knownNB <- Biostrings::translate(DNAStringSet(knownNB_data$Sequence))
#' names(knownNB) <- rownames(knownNB_data)
#' knownNB <- sapply(knownNB, toString)
#' param <- list()
#' param[['NB_Linker1']] <- "GGCCggcggGGCC"
#' param[['NB_Linker2']] <- "GCAGGAGGA"
#' param[['ProteaseSite']] <- "TTAGTCCCAAGA"
#' param[['FC_Linker']] <- "GGCCaaggaggcCGG"
#' param[['knownNB']] <- knownNB
#' param[['nReads']] <- 10000
#' param[['minRelBestHitFreq']] <- 0.8
#' param[['minConsensusScore']] <- 0.9
#' param[['maxMismatch']] <- 1
#' param[['minNanobodyLength']] <- 348
#' param[['minFlycodeLength']] <- 33
#' param[['FCminFreq']] <- 1
#' runNGSAnalysis(file = expFile[1], param)
runNGSAnalysis <- function(file, param){
    param <- .setDefaults(param)
    sampleDir <- NULL
    message(paste0('Read file ', basename(file)))
    myReads <- .getReadsFromFastq(file)
    knownNB <- param[['knownNB']]
    message('... done')
    if (param[['nReads']] > 0 & length(myReads) > param[['nReads']]) {
        myReads <- myReads[seq_len(param[['nReads']])]
    }
    stats <- c(RawReads = length(myReads))
    message(paste0('Filter by ReadLength: ', basename(file)))
    filteredReads <- .filterByReadLengthMedian(myReads, 0.95, 1.05)
    message('... done')
    stats <- c(stats, WidthFilter = length(filteredReads))
    message(paste0('Filter by Pattern: ', basename(file)))
    resultFC <- twoPatternReadFilter(filteredReads,
                                    param[['ProteaseSite']],
                                    param[['FC_Linker']],
                                    maxMismatch = param[['maxMismatch']])
    stats <- c(stats, flankingPatternFC = length(resultFC$reads))
    resultNB <- twoPatternReadFilter(
        resultFC$reads,
        param[['NB_Linker1']],
        param[['NB_Linker2']],
        maxMismatch = param[['maxMismatch']],
        prevPatternPos = resultFC$patternPositions
    )
    stats <- c(stats, flankingPatternNB = length(resultNB$reads))
    mySeq <- list()
    mySeq[['seqFC']] <- subseq(
        resultNB$reads,
        start = resultNB$patternPositions[['leftEnd1']] + 1,
        end = resultNB$patternPositions[['rightStart2']] - 1
    )
    mySeq[['seqNB']] <- subseq(
        resultNB$reads,
        start = resultNB$patternPositions[['leftEnd']] + 1,
        end = resultNB$patternPositions[['rightStart']] - 1
    )
    names(mySeq[['seqFC']]) <- paste('read',
                                    seq(1, length(mySeq[['seqFC']]), 1), sep =
                                        '_')
    names(mySeq[['seqNB']]) <- paste('read',
                                    seq(1, length(mySeq[['seqFC']]), 1), sep =
                                        '_')
    message('... done')
    
    message(paste0('Filter N-Reads: ', basename(file)))
    #Filter reads containing Ns:
    mySeq <- .filterReadsWithNs(mySeq)
    stats <- c(stats, readsWithoutNs = length(mySeq[['seqNB']]))
    message('... done')
    #Filter by width : in frame and flycode >10AS (>32nt)
    
    sampleName <- gsub('.extendedFrags.fastq.gz', '', basename(file))
    nb_width <- sort(table(width(mySeq[['seqNB']])),
                    decreasing = TRUE)[seq_len(10)]
    
    png(paste0('Barplot_NB_Length_', sampleName, '.png'), 600, 400)
    barplot(nb_width[order(names(nb_width))],
            xlab = 'NB length in [nt]',
            ylab = 'Frequency',
            main = 'NB length before Length/InFrame Filtering')
    dev.off()
    
    message(paste0('Filter for in Frame and FlycodeLength: ', basename(file)))
    mySeq <- .filterForInFrame(mySeq,
                            minFlycodeLength = param[['minFlycodeLength']],
                            minNanobodyLength = param[['minNanobodyLength']])
    stats <-
        c(stats, readsWithExpected_FC_NB_Sizes = length(mySeq[['seqNB']]))
    #Translate both to AS:
    mySeqAS <- list()
    mySeqAS[['seqAS_FC']] <- translate(mySeq[['seqFC']])
    mySeqAS[['seqAS_NB']] <- translate(mySeq[['seqNB']])
    message('...done')
    #remove StopCodonSeqs:
    message(paste0('Filter for Stop-Codons: ', basename(file)))
    mySeqAS <- .filterForStopCodons(mySeqAS)
    stats <-
        c(stats, ASSeqWithoutStopCodons = length(mySeqAS[['seqAS_FC']]))
    message('...done')
    #remove ReadsWithWrong Start/End AS:
    message(paste0('Filter for Correct FC Start/End AS Pattern: ',
                    basename(file)))
    mySeqAS <- .filterForStartEndPattern(mySeqAS, myFCPattern = '^GS')
    stats <-
        c(stats, ASSeqWithCorrectFlycodeEnds = length(mySeqAS[['seqAS_FC']]))
    
    bigTable <- data.frame(
        ID = names(mySeqAS[['seqAS_FC']]),
        Flycode = vapply(mySeqAS[['seqAS_FC']],
                        toString, c(FC = '')),
        NanoBody = vapply(mySeqAS[['seqAS_NB']], toString, c(NB = '')),
        stringsAsFactors = FALSE
    )
    
    message('... done')
    #Filter out Flycodes with Freq < minFreq
    message(paste0('Filter for FlyCode-Freq: ', basename(file)))
    mySeqAS_CS <- .filterForMinFlycodeFrequency(
        mySeqAS = mySeqAS,
        minFreq = param[['FCminFreq']],
        file = file,
        knownFC = ''
    )
    
    uniqFCReads <-
        vapply(mySeqAS_CS[['seqAS_FC']], toString, c(FC = ''))
    stats <- c(stats, uniqFC = length(mySeqAS_CS[['seqAS_FC']]))
    message('... done')
    
    ### Export Results:
    message(paste0('Export Result: ', basename(file)))
    
    if(!is.null(sampleDir)){
        if (!file.exists(sampleDir)) {
            dir.create(sampleDir, recursive = TRUE)
        }
        setwd(sampleDir)
    }
     sampleDir <- sub('.fastq.gz', '', basename(file))
    
   
    
    bigTable <- getConsensusNBs(uniqFCReads, bigTable)
    #setwd('..')
    
    bigTable[['FC_SIZE']] <- nchar(bigTable$Flycode)
    bigTable[['NB_SIZE']] <- nchar(bigTable$NanoBody)
    bigTable[['NB_ID']] <- '-'
    bigTable[['CombinationCount']] <- 0
    bigTable[['FC_Count']] <- 0
    bigTable[['NB_Count']] <- 0
    
    for (j in seq_len(nrow(bigTable))) {
        FC_Matches <- (mySeqAS[['seqAS_FC']] == bigTable$Flycode[j])
        bigTable[['FC_Count']][j] <- max(1, sum(FC_Matches))
        NB_Matches <- (mySeqAS[['seqAS_NB']] == bigTable$NanoBody[j])
        bigTable[['NB_Count']][j] <- max(1, sum(NB_Matches))
        CombinationCount <- length(which(FC_Matches & NB_Matches))
        bigTable[['CombinationCount']][j] <- max(1, CombinationCount)
        for (k in seq_len(length(knownNB))) {
            if (bigTable[['NanoBody']][j] == knownNB[[k]])
                bigTable[['NB_ID']][j] <- names(knownNB)[k]
        }
    }
    
    
    bigTableFileName <- basename(sub('.fastq.gz', '_FC2NB.tsv', file))
    write.table(
        bigTable,
        bigTableFileName,
        sep = '\t',
        row.names = FALSE,
        quote = FALSE
    )
    
    
    ##Remove FC where FC count > NB count:
    bigTable <- bigTable[which(bigTable$NB_Count >= bigTable$FC_Count),]
    stats <- c(stats, filter_FC_vs_NB_count = length(unique(bigTable[['Flycode']])))
    
    passMRBHF <- which(bigTable$relBestHitFreq >= param[['minRelBestHitFreq']])
    passCS <- which(bigTable$consensusScore >= param[['minConsensusScore']])
    keepFCs <- rownames(bigTable)[intersect(passMRBHF, passCS)]
    
    stats <- c(stats, FC_consensusFiltering = length(keepFCs))
    
     FC_faFile = file.path(basename(sub('.fastq.gz', '_uniqFC.fasta', file)))
     writeXStringSet(mySeqAS_CS[['seqAS_FC']][keepFCs], file = FC_faFile)
    
    stats <- c(stats, uniqNB = length(unique(bigTable[['NanoBody']])))
    stats <- data.frame(Category = names(stats), stats,
                        stringsAsFactors = FALSE)
    
     statsFileName <- basename(sub('.fastq.gz', '_stats.txt', file))
    
     write.table(
        stats,
        statsFileName,
        sep = '\t',
        quote = FALSE,
        row.names = FALSE
    )
    
    top100_NBFile <- basename(sub('.fastq.gz', '_TopNB.txt', file))
    top100_NB <- .findTopNB(mySeqAS_CS[['seqAS_NB']], n = 100, knownNB)
    
     write.table(
        top100_NB,
        top100_NBFile,
        sep = '\t',
        quote = FALSE,
        row.names = FALSE
    )
    
    #BigTable: Nanobody2Flycode
    uniqNB <- unique(bigTable$NanoBody)
    uniqNB_Summary <-
        data.frame(NB = uniqNB, stringsAsFactors = FALSE)
    uniqNB_Summary[['FlycodeCount']] <- 0
    uniqNB_Summary[['AssociatedFlycodes']] <- ''
    uniqNB_Summary[['NB_Name']] <- ''
    
    for (j in seq_len(length(uniqNB))) {
        flycodes <- bigTable[which(bigTable$NanoBody == uniqNB[j]),][['Flycode']]
        uniqNB_Summary[['FlycodeCount']][j] <- length(flycodes)
        uniqNB_Summary[['AssociatedFlycodes']][j] <- paste(flycodes,
                                                           collapse = ',')
        if (uniqNB[j] %in% knownNB) {
            uniqNB_Summary[['NB_Name']][j] <- 
                names(knownNB)[knownNB == uniqNB[j]]
        }
        uniqNB_Summary <- uniqNB_Summary[order(uniqNB_Summary$NB_Name,
                                            decreasing = TRUE),]
    }
    
     write.table(
        uniqNB_Summary,
        basename(sub('.fastq.gz', '_uniqNB2FC.txt', file)),
        sep = '\t',
        row.names = FALSE,
        quote = FALSE
    )
    
    uniqFC <- unique(bigTable$Flycode)
    uniqFC_Summary <-
        data.frame(FC = uniqFC, stringsAsFactors = FALSE)
    uniqFC_Summary[['NanobodyCount']] <- 0
    uniqFC_Summary[['AssociatedNanobodies']] <- ''
    uniqFC_Summary[['NB_Name']] <- ''
    
    for (j in seq_len(length(uniqFC))) {
        nanobodies <-
            bigTable[which(bigTable$Flycode == uniqFC[j]),][['NanoBody']]
        uniqFC_Summary[['NanobodyCount']][j] <- length(nanobodies)
        concatNBs <- paste(nanobodies, collapse = ',')
        uniqFC_Summary[['AssociatedNanobodies']][j] <- concatNBs
        if (uniqFC_Summary$AssociatedNanobodies[j] %in% knownNB) {
            validNBs <- (knownNB == uniqFC_Summary$AssociatedNanobodies[j])
            uniqFC_Summary[['NB_Name']][j] <- names(knownNB)[validNBs]
        }
    }
    
    orderedNB <- order(uniqFC_Summary$NB_Name, decreasing = TRUE)
    uniqFC_Summary <- uniqFC_Summary[orderedNB, ]
    
    write.table(
        uniqFC_Summary,
        basename(sub('.fastq.gz', '_uniqFC2NB.txt', file)),
        sep = '\t',
        row.names = FALSE,
        quote = FALSE
    )
    
    message('... done')
    class(uniqNB_Summary) <- c(class(uniqNB_Summary), "nanobodyFlycodeLinking")
    return(uniqNB_Summary)
}
