#R

# Christian Panse <cp@fgcz.ethz.ch>
# Lennart Opitz <lopitz@fgcz.ethz.ch>
# 20181128 1547 FGCZ


# Mass spectrometry data are available via ProteomeXchange with identifier
# PXD009301. NGS datasets are available on the European Nucleotide Archive (ENA)
# under accession number PRJEB25673. The NGS and MS data were handled and
# annotated using the B-Fabric information management system28 and are available
# for registered users under project identifiers 1644 and 1875.

NL42_100K <- NestLink:::.getReadsFromFastq("inst/extdata/NL42_100K.fastq.gz")
save(NL42_100K, file="inst/extdata/NestLink_NL42_100K.RData")
