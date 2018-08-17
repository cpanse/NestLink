#' run NGS workflow
#' @param file - sequence file path 
#' @param param - list of input parameters
#' @return character
#' @export runNGSAnalysis
#' @import ShortRead
#' @import Biostrings
runNGSAnalysis = function(file, param){
  # require(ShortRead)
  # require(Biostrings)
  cat(paste0('Read file ', basename(file)))
  myReads = getReadsFromFastq(file)
  cat('...done \n')
  if(param[['nReads']]>0 & length(myReads) > param[['nReads']]){
    myReads = myReads[1:param[['nReads']]]
  }
  stats = c(RawReads=length(myReads))
  cat(paste0('Filter by ReadLength: ', basename(file)))
  filteredReads = filterByReadLengthMedian(myReads,0.95,1.05)
  cat('...done \n')
  stats = c(stats,WidthFilter=length(filteredReads))
  cat(paste0('Filter by Pattern: ', basename(file)))
  resultFC = TwoPatternReadFilter(filteredReads, param[['ProteaseSite']], param[['FC_Linker']], maxMismatch = param[['maxMismatch']])
  stats = c(stats,flankingPatternFC=length(resultFC$reads))
  resultNB = TwoPatternReadFilter(resultFC$reads, param[['NB_Linker1']], param[['NB_Linker2']], maxMismatch = param[['maxMismatch']], 
                                  previousPatternPos=resultFC$patternPositions)
  stats = c(stats,flankingPatternNB=length(resultNB$reads))
  mySeq = list()
  mySeq[['seqFC']] = subseq(resultNB$reads, start=resultNB$patternPositions[['leftEnd1']]+1, end=resultNB$patternPositions[['rightStart2']]-1)
  mySeq[['seqNB']] = subseq(resultNB$reads, start=resultNB$patternPositions[['leftEnd']]+1, end=resultNB$patternPositions[['rightStart']]-1)
  names(mySeq[['seqFC']]) = paste('read',seq(1,length(mySeq[['seqFC']]),1),sep='_')
  names(mySeq[['seqNB']]) = paste('read',seq(1,length(mySeq[['seqFC']]),1),sep='_')
  cat('...done \n')
  
  cat(paste0('Filter N-Reads: ', basename(file)))
  #Filter reads containing Ns:
  mySeq = filterReadsWithNs(mySeq)
  stats = c(stats,readsWithoutNs=length(mySeq[['seqNB']]))
  cat('...done \n')
  #Filter by width : in frame and flycode >10AS (>32nt)
  
  sampleName = gsub('.extendedFrags.fastq.gz','',basename(file))
  nb_width = sort(table(width(mySeq[['seqNB']])),decreasing = T)[1:10]
  png(paste0('Barplot_NB_Length_',sampleName,'.png'), 600, 400)
  barplot(nb_width[order(names(nb_width))],xlab='NB length in [nt]',ylab='Frequency', main='NB length before Length/InFrame Filtering')
  dev.off()
  
  cat(paste0('Filter for in Frame and FlycodeLength: ', basename(file)))
  mySeq = filterForInFrame(mySeq, minFlycodeLength=param[['minFlycodeLength']], minNanobodyLength=param[['minNanobodyLength']])
  stats = c(stats,readsWithExpected_FC_NB_Sizes=length(mySeq[['seqNB']]))
  #Translate both to AS:
  mySeqAS = list()
  mySeqAS[['seqAS_FC']] = translate(mySeq[['seqFC']])
  mySeqAS[['seqAS_NB']] = translate(mySeq[['seqNB']])
  cat('...done \n')
  #remove StopCodonSeqs:
  cat(paste0('Filter for Stop-Codons: ', basename(file)))
  mySeqAS = filterForStopCodons(mySeqAS)
  stats = c(stats,ASSeqWithoutStopCodons=length(mySeqAS[['seqAS_FC']]))
  cat('...done \n')
  #remove ReadsWithWrong Start/End AS:
  cat(paste0('Filter for Correct FC Start/End AS Pattern: ', basename(file)))
  mySeqAS = filterForStartEndPattern(mySeqAS,myFCPattern='^GS')
  stats = c(stats,ASSeqWithCorrectFlycodeEnds=length(mySeqAS[['seqAS_FC']]))
  
  bigTable = data.frame(ID=names(mySeqAS[['seqAS_FC']]),Flycode=sapply(mySeqAS[['seqAS_FC']],toString),
                        NanoBody=sapply(mySeqAS[['seqAS_NB']],toString),stringsAsFactors = F)
  cat('...done \n')
  #Filter out Flycodes with Freq < minFreq
  cat(paste0('Filter for FlyCode-Freq: ', basename(file)))
  mySeqAS_CS = filterForMinFlycodeFrequency(mySeqAS = mySeqAS, minFreq = param[['FCminFreq']], file = file, knownFC = '')
  
  uniqFCReads = sapply(mySeqAS_CS[['seqAS_FC']],toString)
  stats = c(stats,uniqFC=length(mySeqAS_CS[['seqAS_FC']]))
  cat('...done \n')
  
  ###Export Results:
  cat(paste0('Export Result: ', basename(file)))
  setwdNew(gsub('.extended.*','',basename(file)))
  bigTable = getConsensusNBs(uniqFCReads,bigTable)
  setwd('..')
  bigTable[['FC_SIZE']] = nchar(bigTable$Flycode)
  bigTable[['NB_SIZE']] = nchar(bigTable$NanoBody)
  bigTable[['NB_ID']] = '-'
  bigTable[['CombinationCount']] = 0
  bigTable[['FC_Count']] = 0
  bigTable[['NB_Count']] = 0
  
  for(j in 1:nrow(bigTable)){
    bigTable[['FC_Count']][j] = max(1,sum(mySeqAS[['seqAS_FC']]==bigTable$Flycode[j]))
    bigTable[['NB_Count']][j] = max(1,sum(mySeqAS[['seqAS_NB']]==bigTable$NanoBody[j]))
    bigTable[['CombinationCount']][j] = max(1,length(which(mySeqAS[['seqAS_FC']]==bigTable$Flycode[j] & mySeqAS[['seqAS_NB']]==bigTable$NanoBody[j])))
    for(k in 1:length(knownNB)){
      if(bigTable[['NanoBody']][j]==knownNB[[k]])
        bigTable[['NB_ID']][j] = names(knownNB)[k]
    }
  }
  
  
  bigTableFileName = basename(sub('.fastq.gz','_FC2NB.tsv',file))
  write.table(bigTable,bigTableFileName,sep='\t',row.names=F,quote=F)
  
  keepFCs = rownames(bigTable)[bigTable$relBestHitFreq >= param[['minRelBestHitFreq']] & bigTable$consensusScore >= param[['minConsensusScore']]] 
  stats = c(stats,FC_consensusFiltering=length(keepFCs))
  
  FC_faFile = basename(sub('.fastq.gz','_uniqFC.fasta',file))
  writeXStringSet(mySeqAS_CS[['seqAS_FC']][keepFCs],file = FC_faFile) 
  
  stats = c(stats,uniqNB=length(unique(bigTable[['NanoBody']])))
  stats = data.frame(Category=names(stats),stats,stringsAsFactors = F)
  statsFileName = basename(sub('.fastq.gz','_stats.txt',file))
  write.table(stats,statsFileName,sep='\t',quote=F,row.names=F)
  
  top100_NBFile = basename(sub('.fastq.gz','_TopNB.txt',file))
  top100_NB = findTopNB(mySeqAS_CS[['seqAS_NB']],n=100,knownNB)
  write.table(top100_NB,top100_NBFile,sep='\t',quote=F,row.names=F)
  
  #BigTable: Nanobody2Flycode
  uniqNB = unique(bigTable$NanoBody)
  uniqNB_Summary = data.frame(NB = uniqNB, stringsAsFactors = F)
  uniqNB_Summary[['FlycodeCount']] = 0
  uniqNB_Summary[['AssociatedFlycodes']] = ''
  uniqNB_Summary[['NB_Name']] = ''
  for (j in 1:length(uniqNB)){
    flycodes = bigTable[which(bigTable$NanoBody == uniqNB[j]),][['Flycode']]
    uniqNB_Summary[['FlycodeCount']][j] = length(flycodes)
    uniqNB_Summary[['AssociatedFlycodes']][j] = paste(flycodes,collapse=',')
    if(uniqNB[j] %in% knownNB){
      uniqNB_Summary[['NB_Name']][j] = names(knownNB)[which(knownNB == uniqNB[j])]
    }
    uniqNB_Summary = uniqNB_Summary[order(uniqNB_Summary$NB_Name,decreasing =T),]
  }
  write.table(uniqNB_Summary,basename(sub('.fastq.gz','_uniqNB2FC.txt',file)),sep='\t',row.names=F,quote=F)
  
  uniqFC = unique(bigTable$Flycode)
  uniqFC_Summary = data.frame(FC = uniqFC, stringsAsFactors = F)
  uniqFC_Summary[['NanobodyCount']] = 0
  uniqFC_Summary[['AssociatedNanobodies']] = ''
  uniqFC_Summary[['NB_Name']] = ''
  for (j in 1:length(uniqFC)){
    nanobodies = bigTable[which(bigTable$Flycode == uniqFC[j]),][['NanoBody']]
    uniqFC_Summary[['NanobodyCount']][j] = length(nanobodies)
    uniqFC_Summary[['AssociatedNanobodies']][j] = paste(nanobodies,collapse=',')
    if(uniqFC_Summary$AssociatedNanobodies[j] %in% knownNB){
      uniqFC_Summary[['NB_Name']][j] = names(knownNB)[which(knownNB == uniqFC_Summary$AssociatedNanobodies[j])]
    }
  }
  
  uniqFC_Summary = uniqFC_Summary[order(uniqFC_Summary$NB_Name,decreasing = T),]
  write.table(uniqFC_Summary,basename(sub('.fastq.gz','_uniqFC2NB.txt',file)),sep='\t',row.names=F,quote=F)
  cat('...done \n')
  return('success')
}