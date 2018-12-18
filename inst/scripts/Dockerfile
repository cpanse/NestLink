FROM debian:jessie

LABEL version="NestLink 0.99.83"
LABEL description="install NestLink on debian:jessie"

MAINTAINER Christian Panse <cp@fgcz.ethz.ch>

RUN apt-get update && apt-get install subversion gfortran libbz2-dev texlive-extra-utils libpng-dev libblas-dev g++ vim libcurl3 curl -y 
RUN apt-get install `apt-cache depends r-base-dev | grep Depends:|tr "[<>]" " " | awk '{print $NF}' | tr '\n' ' '` -y  
RUN apt-get install  libcurl4-gnutls-dev texlive-extra-utils texlive-latex-recommended texlive-xetex libssl-dev -y

ADD Rdevel.bash /tmp
RUN bash /tmp/Rdevel.bash

RUN echo "BiocManager::install('cpanse/NestLink', version = '3.9', build_vignettes=TRUE)" | /scratch/R-devel//bin/R --no-save
