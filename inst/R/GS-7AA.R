#R

# Christian Panse <cp@fgcz.ethz.ch>
# Pascal Egloff <p.egloff@imm.uzh.ch>

# https://CRAN.R-project.org/package=protViz
# https://github.com/cpanse/NestLink

# [specL](http://dx.doi.org/10.1093/bioinformatics/btv105)
if (packageVersion('protViz') < '0.2.45')
{
  message("protViz package version should be 0.2.45 or higher")
}

library(protViz)
library(gplots)
library(ggplot2)
library(colorspace)

compose_GSx7cTerm <- 
  function(pool=rep("A",7), cTerm=c('WR','WLTVR','WQEGGR','WQSR','WLR')){ 
  paste("GS", 
        paste(pool[sample(length(pool), 7)], collapse=''), 
        cTerm[sample(length(cTerm), 1)], 
        sep='') 
  }


## Compose a GPYYXXXXXXYYR peptide
compose_GPx10R <- function(aa_pool1, aa_pool2){ 
  paste("GP", paste(aa_pool1[sample(length(aa_pool1), 2)], collapse=''), 
        paste(aa_pool2[sample(length(aa_pool2), 6)], collapse=''),
        paste(aa_pool1[sample(length(aa_pool1), 2)], collapse=''), "R",
        sep='') 
  }

hydrophobicity <- function(x){ ssrc(x) }

.in_silico_LCMS_map <- function(x, bins=10, ...){
  hyd <- hydrophobicity(x)
  pim <- unlist(lapply(x, function(x){parentIonMass(x)}))

  df <- data.frame(hyd=hyd, pim=pim)
  
  p <- ggplot(df, aes(hyd, pim)) + 
    stat_bin2d(bins=80) +  
    labs(title = "in-silico LCMS map", 
         subtitle = paste(deparse(substitute(x)), "| sample size =", length(x))) + 
    labs(x = "hydrophobicity value [as computed by SSRC; dimensionless quantity]", y = "parent ion mass [in Dalton]") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    scale_fill_gradientn(colours = rev(colorspace::sequential_hcl(8)))
  p
}

#######################################################################################################
## Input
(sample.size <- 10000)
aa_pool_x7 <- c(rep('A', 18), rep('S', 6), rep('T', 12), rep('N', 1), rep('Q', 1), rep('D', 11), 
              rep('E', 11), rep('V', 12), rep('L', 2), rep('F', 1), rep('Y', 4), rep('W', 1), 
              rep('G', 8), rep('P', 12))

mycTerm <- c('WR','WLTVR','WQEGGR','WQSR','WLR')

## Compose a GSXXXXXXX(WR|WLTVR|WQGGER|WQSR|WLR) peptide
set.seed(2)
peptides.GSx7cTerm <- replicate(sample.size, compose_GSx7cTerm(pool=aa_pool_x7, cTerm=mycTerm))


## Some Sanity Checks
table(aa_pool_x7)

stopifnot(length(aa_pool_x7) == 100)

stopifnot(length(peptides.GSx7cTerm[grepl("^GS[ASTNQDEFVLYWGP]{7}(WR|WLTVR|WQEGGR|WLR|WQSR)$", peptides.GSx7cTerm)]) 
          == sample.size)



## MAIN
pdf("~/NestLink.pdf", 12,10)

print(.in_silico_LCMS_map(peptides.GSx7cTerm, main='GSx7cTerm'))

dev.off()
