# NestLink

Engineered Peptide Barcodes for In-Depth Analyses of Binding Protein Ensembles


## 1. System requirements


### Software dependencies
- install R (> 3.4.0)

- pandoc (for exporting markdown files to other formats).

- install https://CRAN.R-project.org/package=devtools

- install required R packages

**Whatsoever for a clean install, I recommend using the current R version 3.5 and the current Bioconductor 3.7.**


### Versions the software has been tested on

|platform|platform version|R version|note|
| :------- |:---------------| -------:|:------- |
|Linux     | Debian 10 ([buster](https://www.debian.org/releases/testing/releasenotes)) | 3.5.0 , Bioconductor version 3.7| CP |
|Microsoft |Server 2012 R2 x64| R (3.5.0), Bioconductor version 3.7||
| macOS High| Sierra 10.13.4| R (3.4.2)||


## 2. Installation guide

run an R session and execute the following R code snippet

```{r}
library(devtools)
install_github('cpanse/NestLink', build_vignettes = TRUE)
```

alternative code snippets in case the line above fails.
```{r}
https://git.bioconductor.org/packages/specL
#
source("https://bioconductor.org/biocLite.R")
biocLite("specL")
biocLite("BiocStyle")
#
install_git('https://git.bioconductor.org/packages/specL')

```

**Typical install time** - 
Expect an hour to get all the R packages (BioConductor) running. 
If all dependencies are installed 
`install_github('cpanse/NestLink', build_vignettes = TRUE)` requires 1m47.014s 
on a Intel(R) Core(TM) i5-2500 CPU @ 3.30GHz running Debian Linux.

## 3. Demonstration / Documentation

Instructions to run on data and expected output is described in the package's 
vignettes.

```{r}
browseVignettes('NestLink')
```

Expected run time for the vignette build is 1 minute on a today's desktop computer.

## 4. Instructions for use

read the vignettes.

```{r}
browseVignettes('NestLink')
```

## References 

- [project p1875 at the Functional Genomics Center Zurich](https://fgcz-bfabric.uzh.ch/bfabric/userlab/show-project.html?id=1875)

- https://www.biorxiv.org/content/early/2018/03/23/287813


