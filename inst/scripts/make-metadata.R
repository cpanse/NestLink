#R

### =========================================================================
### NestLink metadata 
### -------------------------------------------------------------------------
###

# Christian Panse <cp@fgcz.ethz.ch>
# Lennart Opitz <lopitz@fgcz.ethz.ch>
# 20181128 1347 FGCZ


meta <- data.frame(
    Title = c(paste0("Sample data NGS"),
              paste0("Flycode tryptic digested"),
              paste0("Nanobody tryptic digested"),
              paste0("FASTA as ground-truth for unit testing"),
              paste0("Known nanobodies"),
              paste0("Quantitaive results for SMEG and COLI")
              ),
    Description = c(paste0("sample data NGS in FASTQ"),
                    paste0("Flycode tryptic digested AA sequences with ESP_Prediction"),
                    paste0("Nanobody tryptic digested AA sequences with ESP_Prediction"),
                    paste0("FASTA data as ground-truth for unit testing"),
                    paste0("known nanobodies in NC sequences"),
                    paste0("mass spectrometry quantitative results of nanobodies expressed in SMEG and COLI species.")
            ),
    BiocVersion = c("3.8", "3.9"),
    Genome = NA, 
    SourceType = c('fastq.gz', 'txt', 'txt', 'RData', 'txt', 'csv'),
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
    DataProvider = "FGCZ",
    Maintainer = c("Lennart Opitz <lopitz@fgcz.ethz.ch>",
                   "Christian Panse <cp@fgcz.ethz.ch>",
                   "Christian Panse <cp@fgcz.ethz.ch>",
                   "Lennart Opitz <lopitz@fgcz.ethz.ch>",
                   "Lennart Opitz <lopitz@fgcz.ethz.ch>",
                   "Christian Panse <cp@fgcz.ethz.ch>"),
    RDataClass = c("DNAStringSet", "data.frame", "data.frame", "data.frame","data.frame","data.frame"),
    DispatchClass = NA,
    RDataPath = c("NestLink/NL42_100K.fastq.gz",
                  "NestLink/FC.tryptic",
                  "NestLink/NB.tryptic",
                  "NestLink/nanobodyFlycodeLinkage.RData",
                  "NestLink/knownNB.txt",
                  "NestLink/PGexport2_normalizedAgainstSBstandards_Peptides.csv"
                  ),
    Tags = "",
    Notes = NA
)

write.csv(meta, file="inst/extdata/metadata.csv", row.names=FALSE)
    