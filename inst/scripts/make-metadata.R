#R

### =========================================================================
### NestLink metadata 
### -------------------------------------------------------------------------
###

# Christian Panse <cp@fgcz.ethz.ch>
# Lennart Opitz <lopitz@fgcz.ethz.ch>
# 20181128 1347 FGCZ


meta <- data.frame(
    Title = c(paste0("Sample NGS NB FC linkage data"),
              paste0("Flycode tryptic digested"),
              paste0("Nanobody tryptic digested"),
              paste0("FASTA as ground-truth for unit testing"),
              paste0("Known nanobodies"),
              paste0("Quantitaive results for SMEG and COLI")
              ),
    Description = c(paste0("Sample NGS demonstratig the linkage between nanobodies (NB) and flycodes (FC). data in FASTQ"),
                    paste0("Flycode tryptic digested amino acid sequences with ESP_Prediction score."),
                    paste0("Nanobody tryptic digested amino acid sequences with ESP_Prediction score."),
                    paste0("FASTA data as ground-truth for unit testing"),
                    paste0("Known nanobodies as nucleic acid sequences"),
                    paste0("Mass spectrometry based quantitative results of nanobodies expressed in SMEG and COLI species.")
            ),
    BiocVersion = rep("3.9", 6),
    Genome = NA, 
    SourceType = c('FASTQ', 'TXT', 'TXT', 'RData', 'TXT', 'CSV'),
    SourceUrl = c("https://fgcz-bfabric.uzh.ch/bfabric/userlab/show-project.html?id=1644",
                  "https://fgcz-bfabric.uzh.ch/bfabric/userlab/show-project.html?id=1875",
                  "https://fgcz-bfabric.uzh.ch/bfabric/userlab/show-project.html?id=1875",
                  "https://fgcz-bfabric.uzh.ch/bfabric/userlab/show-project.html?id=1644",
                  "https://fgcz-bfabric.uzh.ch/bfabric/userlab/show-project.html?id=1644",
                  "https://fgcz-bfabric.uzh.ch/bfabric/userlab/show-project.html?id=1875"),
    SourceVersion = "Nov 28 2018",
    Species = NA,
    TaxonomyId = NA,
    Coordinate_1_based = NA,
    DataProvider = "Functional Genomics Center Zurich",
    Maintainer = c("Lennart Opitz <lopitz@fgcz.ethz.ch>",
                   "Christian Panse <cp@fgcz.ethz.ch>",
                   "Christian Panse <cp@fgcz.ethz.ch>",
                   "Lennart Opitz <lopitz@fgcz.ethz.ch>",
                   "Lennart Opitz <lopitz@fgcz.ethz.ch>",
                   "Christian Panse <cp@fgcz.ethz.ch>"),
    RDataClass = c("DNAStringSet", "data.frame", "data.frame", "data.frame","data.frame","data.frame"),
    DispatchClass = rep("FilePath", 6),
    RDataPath = c("NestLink/NL42_100K.fastq.gz",
                  "NestLink/FC.tryptic",
                  "NestLink/NB.tryptic",
                  "NestLink/nanobodyFlycodeLinkage.RData",
                  "NestLink/knownNB.txt",
                  "NestLink/PGexport2_normalizedAgainstSBstandards_Peptides.csv"
                  ),
    Tags = "",
    Notes = c(NA, NA, NA, NA, NA, "bfabric WU158716, WU158717"),
    MD5 = c("4a13c5c61a5b29f4fd8830c1c15419b6",
            "f6faa7458350ce1805bec30e9ffdeaae",
            "db85a806c5151113536b710d566d9cf3",
            "57b2756fb0ebcf73d4036846580cb5b2",
            "003bf82c58f0a96a2bd945d171dc907c",
            "0ca525d0a65d4938f0cbc785b7e0d2d3")
)

write.csv(meta, file="inst/extdata/metadata.csv", row.names=FALSE)
    
