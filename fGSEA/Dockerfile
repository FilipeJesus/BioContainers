FROM r-base:3.5.0

# File Author / Maintainer 
MAINTAINER Filipe Jesus

# Install bioconductor 
RUN R -e 'source("http://bioconductor.org/biocLite.R"); biocLite("fgsea", ask= FALSE, version = "3.8")'

# Cleanup 
RUN apt-get clean
