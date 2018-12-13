#R

### =========================================================================
### NestLink metadata 
### -------------------------------------------------------------------------
###

# Christian Panse <cp@fgcz.ethz.ch>
# Lennart Opitz <lopitz@fgcz.ethz.ch>
# 20181128 1347 FGCZ
# 20181212 0852 FGCZ

n <- 8

meta <- data.frame(
    Title = c(paste0("Sample NGS NB FC linkage data"),
              paste0("Flycode tryptic digested"),
              paste0("Nanobody tryptic digested"),
              paste0("FASTA as ground-truth for unit testing"),
              paste0("Known nanobodies"),
              paste0("Quantitaive results for SMEG and COLI"),
              paste0("F255744 Mascot Search result"),
              paste0("WU160118 Mascot Search results")
              ),
    Description = c(paste0("Sample NGS demonstratig the linkage between nanobodies (NB) and flycodes (FC). data in FASTQ"),
                    paste0("Flycode tryptic digested amino acid sequences with ESP_Prediction score."),
                    paste0("Nanobody tryptic digested amino acid sequences with ESP_Prediction score."),
                    paste0("FASTA data as ground-truth for unit testing"),
                    paste0("Known nanobodies as nucleic acid sequences"),
                    paste0("Mass spectrometry based quantitative results of nanobodies expressed in SMEG and COLI species."),
                    paste0("F255744 Mascot Search result"),
                    paste0("WU160118 Mascot Search results")
            ),
    BiocVersion = rep("3.9", n),
    Genome = rep(NA, n), 
    SourceType = c('FASTQ', 'TXT', 'TXT', 'RData', 'TXT', 'CSV', 'TXT', 'TXT'),
    SourceUrl = c("https://fgcz-bfabric.uzh.ch/bfabric/userlab/show-project.html?id=1644",
                  "https://fgcz-bfabric.uzh.ch/bfabric/userlab/show-project.html?id=1875",
                  "https://fgcz-bfabric.uzh.ch/bfabric/userlab/show-project.html?id=1875",
                  "https://fgcz-bfabric.uzh.ch/bfabric/userlab/show-project.html?id=1644",
                  "https://fgcz-bfabric.uzh.ch/bfabric/userlab/show-project.html?id=1644",
                  "https://fgcz-bfabric.uzh.ch/bfabric/userlab/show-project.html?id=1875",
                  "https://fgcz-bfabric.uzh.ch/bfabric/userlab/show-resource.html?id=409912", 
                  "https://fgcz-bfabric.uzh.ch/bfabric/userlab/show-workunit.html?id=160118"),
    SourceVersion = c(rep("Nov 28 2018", 6), rep("Dec 13 2018",2)),
    Species = rep(NA, n),
    TaxonomyId = rep(NA, n),
    Coordinate_1_based = rep(NA, n),
    DataProvider = "Functional Genomics Center Zurich",
    Maintainer = c("Lennart Opitz <lopitz@fgcz.ethz.ch>",
                   "Christian Panse <cp@fgcz.ethz.ch>",
                   "Christian Panse <cp@fgcz.ethz.ch>",
                   "Lennart Opitz <lopitz@fgcz.ethz.ch>",
                   "Lennart Opitz <lopitz@fgcz.ethz.ch>",
                   "Christian Panse <cp@fgcz.ethz.ch>",
                   "Christian Panse <cp@fgcz.ethz.ch>",
                   "Christian Panse <cp@fgcz.ethz.ch>"),
    RDataClass = c("DNAStringSet", rep("data.frame", 7)),
    DispatchClass = rep("FilePath", n),
    RDataPath = c("NestLink/NL42_100K.fastq.gz",
                  "NestLink/FC.tryptic",
                  "NestLink/NB.tryptic",
                  "NestLink/nanobodyFlycodeLinkage.RData",
                  "NestLink/knownNB.txt",
                  "NestLink/PGexport2_normalizedAgainstSBstandards_Peptides.csv",
                  "NestLink/F255744.RData",
                  "NestLink/WU160118.RData"
                  ),
    Tags = "",
    Notes = paste("Markus Seeger <m.seeger@imm.uzh.ch>, Pascal Egloff <p.egloff@imm.uzh.ch>", 
                  c("", "", "", "", "", " bfabric WU158716, WU158717", "" , " bfabric WU160118"), sep=';'),
    MD5 = c("4a13c5c61a5b29f4fd8830c1c15419b6",
            "f6faa7458350ce1805bec30e9ffdeaae",
            "db85a806c5151113536b710d566d9cf3",
            "57b2756fb0ebcf73d4036846580cb5b2",
            "003bf82c58f0a96a2bd945d171dc907c",
            "0ca525d0a65d4938f0cbc785b7e0d2d3",
            "d5e4d13e9ecba4231d1808c6bb0bb454",
            "a17f4505e322d440bc0e9edf8e5277bb"
            )
)

write.csv(meta, file="inst/extdata/metadata.csv", row.names=FALSE)
    
