FROM r-base:3.5.0

# File Author / Maintainer
MAINTAINER Filipe Jesus

# Install bioconductor
RUN apt-get update
RUN apt-get install -y aptitude
RUN apt-get install -y libcurl4-openssl-dev
RUN apt-get install -y libxml2-dev
RUN apt-get install -y libssl-dev
RUN R -e 'source("http://bioconductor.org/biocLite.R"); biocLite(c("GenomicRanges","GenomeInfoDb","Rsamtools","GenomicAlignments","GenomicFeatures","dplyr","rtracklayer", "Rsubread"), ask= FALSE, version = "3.8")'

# Cleanup
RUN apt-get clean
