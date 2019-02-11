FROM r-base:3.5.0

# File Author / Maintainer 
MAINTAINER Filipe Jesus

# Install bioconductor 
RUN R -e 'source("http://bioconductor.org/biocLite.R"); biocLite("fgsea")'

# Cleanup 
RUN apt-get clean
