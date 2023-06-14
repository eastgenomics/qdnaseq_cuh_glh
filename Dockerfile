## R packages

FROM rocker/r-ver:4.3.0

RUN mkdir /home/software
# install curl and libcurl
RUN apt update
RUN apt-get install -y curl
RUN apt-get -y install libcurl4-openssl-dev
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils
RUN apt-get -y update && apt-get install -y \
   default-jdk  r-cran-rjava  r-cran-nloptr libssh2-1-dev  libz-dev

# install R packages

RUN Rscript -e "install.packages('BiocManager', dependencies = TRUE)"
RUN Rscript -e "install.packages('devtools')"
RUN Rscript -e "install.packages('tidyverse')"
RUN Rscript -e "install.packages('Matrix', dependencies = TRUE)"
RUN Rscript -e "BiocManager::install(c('QDNAseq', 'MatrixGenerics', 'IRanges', 'S4Vectors', 'GenomicRanges', 'cn.mops', 'ACE', 'BSgenome.Hsapiens.UCSC.hg38'))"
RUN Rscript -e 'devtools::install_github(repo = "asntech/QDNAseq.hg38")'

WORKDIR /home/software
