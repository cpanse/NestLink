#!/bin/bash
#
# Installs R-devel  
#
# Christian Panse <cp@fgcz.ethz.ch>
# $HeadURL: svn+ssh://cp@fgcz-148.uzh.ch/home/cp/__SUBVERSION_REPOSITORY__/__projects/2016/20160701__daily_Rdevel_build.bash $
# $Id: 20160701__daily_Rdevel_build.bash 1126 2018-01-31 05:05:26Z cp $
# $Date: 2018-01-31 06:05:26 +0100 (Wed, 31 Jan 2018) $

set -x
set -e
set -o pipefail


RDEVELROOT=/scratch/R-devel/

mkdir -p $RDEVELROOT \
  && cd $RDEVELROOT \
  && svn co https://svn.r-project.org/R/trunk $RDEVELROOT/source --trust-server-cert \
  && cd $RDEVELROOT \
  && cd source \
  && ./configure --prefix=$RDEVELROOT --without-recommended-packages --enable-R-shlib --with-x=no  \
  && time nice -19 make -j 32 && make install 

$RDEVELROOT/bin/R --no-save <<EOF
#R
install.packages(c('remotes', 'BiocManager'), repos='http://stat.ethz.ch/CRAN')
EOF

echo "BiocManager::install('cpanse/NestLink', version = '3.9')" \
  | $RDEVELROOT/bin/R --no-save 
