# NestLink

Engineered Peptide Barcodes for In-Depth Analyses of Binding Protein Ensembles


## 1. System requirements


### Software dependencies
- install R (> 3.4.0)

- install https://CRAN.R-project.org/package=devtools

- install required R packages

### Versions the software has been tested on

- Debian Linux, R (3.5)

- Windows 10, R (3.4.4)

**Whatsoever for a clean install, I recommend using the current R version 3.5 and the current Bioconductor 3.7.**


## 2. Installation guide

run a R session and execute the following R code snippet

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite("specL")
biocLite("BiocStyle")

library(devtools)
install_github('cpanse/NestLink', build_vignettes = TRUE)
```

**Typical install time**
Expect an hour to get all the R packages running.

## 3. Demo 

Instructions to run on data and expected output is described in the package's 
vignettes.

```{r}
browseVignettes('NestLink')
```

Expected run time for the vignette build is 10 minutes on a todays desktop computer.

## 4. Instructions for use

can be found in the vignettes.
```{r}
browseVignettes('NestLink')
```

## References 

- [FGCZ project p1875  NestLink](https://fgcz-bfabric.uzh.ch/bfabric/userlab/show-project.html?id=1875)

- https://www.biorxiv.org/content/early/2018/03/23/287813


