FROM debian:jessie

LABEL version="NestLink 0.99.83"
LABEL description="install NestLink on debian:jessie"

MAINTAINER Christian Panse <cp@fgcz.ethz.ch>

ADD Rdevel.bash /tmp
RUN apt-get update && apt-get install subversion gfortran libbz2-dev texlive-extra-utils libpng-dev libblas-dev g++ vim -y 
RUN do apt-get install `apt-cache depends r-base-dev|awk '{print $NF}' | tr "\n" " "` -y  
RUN bash /tmp/Rdevel.bash
