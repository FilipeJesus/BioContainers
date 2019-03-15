FROM machalen/dexseq

# File Author / Maintainer 
MAINTAINER Filipe Jesus

# Install bioconductor 
RUN R -e "install.packages('dplyr')"

# Cleanup 
RUN apt-get clean
