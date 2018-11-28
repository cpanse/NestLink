#R

# Christian Panse <cp@fgcz.ethz.ch>
# Lennart Opitz <lopitz@fgcz.ethz.ch>
# 20181128 1547 FGCZ

NL42_100K <- NestLink:::.getReadsFromFastq("inst/extdata/NL42_100K.fastq.gz")


save(NL42_100K, 
    file="NestLink_NL42_100K.Rda")
