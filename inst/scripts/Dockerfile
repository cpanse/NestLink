FROM debian:jessie

LABEL version="NestLink 0.99.83"
LABEL description="install NestLink on debian:jessie"

MAINTAINER Christian Panse <cp@fgcz.ethz.ch>

RUN apt-get update \
  && apt-get install subversion gfortran libbz2-dev texlive-extra-utils libpng-dev libblas-dev g++ vim libcurl3 curl -y 
RUN apt-get install `apt-cache depends r-base-dev | grep Depends:|tr "[<>]" " " | awk '{print $NF}' | tr '\n' ' '` -y  
RUN apt-get install libcurl4-gnutls-dev texlive-extra-utils texlive-latex-recommended texlive-xetex libssl-dev -y

RUN RDEVELROOT=/scratch/R-devel/ \
  && mkdir -p $RDEVELROOT \
  && cd $RDEVELROOT \
  && svn --quiet checkout https://svn.r-project.org/R/trunk $RDEVELROOT/source --trust-server-cert \
  && cd $RDEVELROOT \
  && cd source \
  && ./configure --prefix=$RDEVELROOT --without-recommended-packages --enable-R-shlib --with-x=no  \
  && make -j 32 \
  && make install \
  && echo "install.packages(c('remotes', 'BiocManager'), repos='http://stat.ethz.ch/CRAN')" \
  | $RDEVELROOT/bin/R --no-save 

RUN echo "BiocManager::install('cpanse/NestLink', version = '3.9', build_vignettes=FALSE)" \
  | /scratch/R-devel/bin/R --no-save
