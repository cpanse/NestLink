# NestLink

Engineered Peptide Barcodes for In-Depth Analyses of Binding Protein Ensembles

## 1. System requirements

### Software dependencies

- install R (>= 3.6)

- install Bioconductor (>=3.9)

## 2. Installation guide

run an R session and execute the following R code snippet

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install("cpanse/NestLink", version = "3.9")  
```

or using docker

```
docker pull cpanse/nestlink \
  && docker run -a stdin -a stdout -i -t cpanse/nestlink /scratch/R-devel/bin/R
```


**Typical install time** - 
based on the [Dockerfile](inst/scripts/Dockerfile) the install snippet above 
took 19m46.464s on a linux server (RAID6, Intel(R) Xeon(R) CPU E5-2698 v3 @ 2.30GHz) and one hour on dockerhub.

As an alternative you can also consider the latest [release](https://github.com/cpanse/NestLink/releases).

### Versions the software has been tested on

|platform|NestLink version|platform version|R version|note|
| :------- |:---------------|:---------------| :-------|:------- |
|Linux     | 0.99.51 | Debian 10 ([buster](https://www.debian.org/releases/testing/releasenotes)) | R 3.5.1, Bioconductor version 3.8| CP |
|Microsoft | 0.99.51 | Server 2012 R2 x64| R 3.5.0, Bioconductor version 3.7||
| macOS High| 0.99.51 | Sierra 10.13.4| R 3.4.2||

## 3. Demonstration / Documentation

Instructions to run on data and expected output is described in the package's 
vignettes.

```{r}
browseVignettes('NestLink')
```

please study the vignettes in the following order:

0. Derive Peptide FlyCodes by Conducting Random Experiment  
1. NGS filtering workflow to get high quality FlyCode and Nanobody sequences  
2. FASTA p1875 db10 - ESP / SSRC prediction - Summary 
3. Compare Predicted and Measured FlyCodes (F255744).  
4. Control experiment to assess robustness of protein detection via flycodes 


Expected run time for the vignette build is less than 5 minutes on a today's desktop computer.

## 4. Instructions for use

read the vignettes.

```{r}
browseVignettes('NestLink')
```

## References 

- [project p1875 at the Functional Genomics Center Zurich](https://fgcz-bfabric.uzh.ch/bfabric/userlab/show-project.html?id=1875)

- [bioRxiv 2018/03/23/287813](https://www.biorxiv.org/content/early/2018/03/23/287813)


