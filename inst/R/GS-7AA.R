#R

# Christian Panse <cp@fgcz.ethz.ch>
# Bernd Roschitzki 
# Pascal Egloff <p.egloff@imm.uzh.ch>
# 2015 - 2018

# https://CRAN.R-project.org/package=protViz
# https://github.com/cpanse/NestLink

stopifnot(packageVersion('protViz') >= '0.2.45')
stopifnot(packageVersion('NestLink') >= '0.99.14')


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

hydrophobicity <- 
  function(x){ unlist(lapply(x, ssrc)) }

.in_silico_LCMS_map <- 
  function(x, bins=10, alpha = 0.5, size = 1.0, ...){
  hyd <- hydrophobicity(x)
  pim <- unlist(lapply(x, function(x){parentIonMass(x)}))

  df <- data.frame(hyd=hyd, pim=pim)
  
  p <- ggplot(df, aes(hyd, pim)) + 
    geom_point(color=rgb(1, 0.6, 0.6, alpha = alpha), size = size, pch = 16) +
    labs(title = "in-silico LC-MS map", 
         subtitle = paste(deparse(substitute(x)), "| sample size =", length(x), "| size =", size)) + 
    labs(x = "hydrophobicity value [as computed by SSRC; dimensionless quantity]",
         y = "peptide mass [in Dalton]") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    scale_fill_gradientn(colours = diverge_hsv(32, h = 0, s = c(1, 0.4), v = c(0.8, 1.0), alpha = 0.5)[17:(17+8)])
  p
}

#######################################################################################################
## Input
(sample.size <- 10000)
aa_pool_x7 <- c(rep('A', 18), rep('S', 6), rep('T', 12), rep('N', 1), rep('Q', 1), rep('D', 11), 
              rep('E', 11), rep('V', 12), rep('L', 2), rep('F', 1), rep('Y', 4), rep('W', 1), 
              rep('G', 8), rep('P', 12))


## Compose a GSXXXXXXX(WR|WLTVR|WQGGER|WQSR|WLR) peptide
set.seed(2)
FC.GSx7cTerm <- replicate(sample.size, compose_GSx7cTerm(pool=aa_pool_x7))

## Some Sanity Checks
table(aa_pool_x7)

stopifnot(length(aa_pool_x7) == 100)

stopifnot(length(FC.GSx7cTerm[grepl("^GS[ASTNQDEFVLYWGP]{7}(WR|WLTVR|WQEGGR|WLR|WQSR)$", FC.GSx7cTerm)]) 
          == sample.size)




## MAIN
pdf("~/NestLink.pdf", 6, 0.6 * 6)
op<-par(mfrow=c(2,1))
lapply (seq(0.2, 4, by=0.2), function(size){print(.in_silico_LCMS_map(FC.GSx7cTerm, alpha = 0.5, size = size))})


ESP <- rbind(.getFC(), .getNB())
p <- ggplot(ESP, aes(x = ESP_Prediction, fill = cond)) +
  geom_histogram(bins = 50, alpha = 0.4, position="identity") +
  labs(x = "detectability in LC-MS (ESP prediction)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(p)

ESP <- rbind(.getFC(), NB.unambiguous(.getNB()))
p <- ggplot(ESP, aes(x = ESP_Prediction, fill = cond)) +
  geom_histogram(bins = 50, alpha = 0.4, position="identity") +
  labs(x = "detectability in LC-MS (ESP prediction)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(p)

ESP <- rbind(.getFC(), NB.unique(.getNB()))
p <- ggplot(ESP, aes(x = ESP_Prediction, fill = cond)) +
  geom_histogram(bins = 50, alpha = 0.4, position="identity") +
  labs(x = "detectability in LC-MS (ESP prediction)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(p)

  
ESP <- rbind(.getNB(), NB.unambiguous(.getNB()), NB.unique(.getNB()))

p <- ggplot(ESP, aes(x = ESP_Prediction, fill = cond)) +
  geom_histogram(bins = 50, alpha = 0.4, position="identity") +
  labs(x = "detectability in LC-MS (ESP prediction)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(p)

table(ESP$cond)

dev.off()


