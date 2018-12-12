#R

context("Test NGS")

test_that("NGS ground truth", {
  
  library(ExperimentHub)
  
  eh <- ExperimentHub()
  
  f <- query(eh, c("NestLink", "nanobodyFlycodeLinkage.RData"))[[1]]
  
  # f <- system.file("extdata/nanobodyFlycodeLinkage.RData", package="NestLink")
  expect_true(file.exists(f))
  load(f)
  expect_true(NestLink:::is.nanobodyFlycodeLinking(nanobodyFlycodeLinkage.sample))
  
})


test_that("NGS pipeline", {
  library(ExperimentHub)
  
  eh <- ExperimentHub()
  
  f <- query(eh, c("NestLink", "nanobodyFlycodeLinkage.RData"))[[1]]
  
   # f <- system.file("extdata/nanobodyFlycodeLinkage.RData", package="NestLink")
    expect_true(file.exists(f))
    load(f)
    expect_true(NestLink:::is.nanobodyFlycodeLinking(nanobodyFlycodeLinkage.sample))
    
    # expFile <-system.file("extdata/NL42_100K.fastq.gz", package="NestLink")
    expFile <- query(eh, c("NestLink", "NL42_100K.fastq.gz"))[[1]]
    expect_true(file.exists(expFile))
    scratchFolder <- tempdir()
    setwd(scratchFolder)
    
    
    # knownNB_File <- system.file("extdata/knownNB.txt", package="NestLink")
    knownNB_File <- query(eh, c("NestLink", "knownNB.txt"))[[1]]
    knownNB_data <- read.table(knownNB_File,
                               sep='\t',
                               header = TRUE,
                               row.names = 1,
                               stringsAsFactors = FALSE)
    
    knownNB <- Biostrings::translate(DNAStringSet(knownNB_data$Sequence))
    names(knownNB) <- rownames(knownNB_data)
    knownNB <- sapply(knownNB, toString)
    
    param <- list()
    param[['NB_Linker1']] <- "GGCCggcggGGCC"
    param[['NB_Linker2']] <- "GCAGGAGGA"
    param[['ProteaseSite']] <- "TTAGTCCCAAGA"
    param[['FC_Linker']] <- "GGCCaaggaggcCGG"
    param[['knownNB']] <- knownNB
    param[['nReads']] <- 100
    param[['minRelBestHitFreq']] <- 0.8 
    param[['minConsensusScore']] <- 0.9
    param[['maxMismatch']] <- 1
    param[['minNanobodyLength']] <- 348
    param[['minFlycodeLength']] <- 33
    param[['FCminFreq']] <- 1
    
    NB2FC <- runNGSAnalysis(file = expFile[1], param)
    
    
    expect_true(sum((nanobodyFlycodeLinking.summary(NB2FC) == 
                         nanobodyFlycodeLinking.summary(nanobodyFlycodeLinkage.sample))) == 3)
})


