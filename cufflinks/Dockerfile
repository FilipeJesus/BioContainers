#Build the image based on Ubuntu
FROM python:2.7.15

#Maintainer and author
MAINTAINER Filipe Jesus <filipeje@hotmail.co.uk>
#Insperation from Magdalena Arnal <marnal@imim.es> docker machalen/cufflinksdocker

ADD cufflinks/tablemaker2.1.1.Linux_x86_64.tar.gz /

#Update and install packages wget, unzip
RUN apt-get update -y && apt-get install -y \
    wget git unzip bzip2 g++ make zlib1g-dev ncurses-dev default-jdk default-jre libncurses5-dev \
    libbz2-dev liblzma-dev

#Set wokingDir in /bin
WORKDIR /bin


#Download Cufflinks
RUN wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz
#Unzip Cufflinks
RUN tar zxvf cufflinks-2.2.1.Linux_x86_64.tar.gz
#Clean
RUN rm cufflinks-2.2.1.Linux_x86_64.tar.gz

RUN pip install numpy
RUN pip install pandas
RUN pip install scipy 
RUN pip install statsmodels

#Download samtools from GitHub
RUN wget http://github.com/samtools/samtools/releases/download/1.5/samtools-1.5.tar.bz2
RUN tar --bzip2 -xf samtools-1.5.tar.bz2
WORKDIR /bin/samtools-1.5
RUN ./configure
RUN make
RUN rm /bin/samtools-1.5.tar.bz2

#Add Cufflinks and samtools to the path variable
ENV PATH $PATH:/bin/cufflinks-2.2.1.Linux_x86_64
ENV PATH $PATH:/bin/samtools-1.5

#Set the default Working Directory
WORKDIR /

# Install tablemaker
RUN apt-get update
RUN cp -r tablemaker-2.1.1.Linux_x86_64/ /bin/
RUN rm -rf tablemaker-2.1.1.Linux_x86_64/

ENV PATH $PATH:/bin/tablemaker-2.1.1.Linux_x86_64
