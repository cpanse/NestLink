# NestLink

Engineered Peptide Barcodes for In-Depth Analyses of Binding Protein Ensembles


## 1. System requirements

- install R (> 3.4.0)

- install https://CRAN.R-project.org/package=devtools

- install required R packages



## 2. Installation guide

run a R session and execute the following R code snippet
```{r}
source("https://bioconductor.org/biocLite.R")
biocLite("specL")

library(devtools)
install_github('cpanse/NestLink', build_vignettes = TRUE)
```

## 3. Demo 

- [vignettes/DeriveFlyCodes.Rmd](vignettes/DeriveFlyCodes.Rmd)

```{r}
browseVignettes('NestLink')
```


## 4. Instructions for use


## References 

- [FGCZ project p1875  NestLink](https://fgcz-bfabric.uzh.ch/bfabric/userlab/show-project.html?id=1875)

- https://www.biorxiv.org/content/early/2018/03/23/287813


