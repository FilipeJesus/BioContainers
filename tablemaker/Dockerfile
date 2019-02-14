#Build the image based on machalen's cufflinks can be found on Dockerhub
FROM machalen/cufflinksdocker:latest

# File Author / Maintainer
MAINTAINER Filipe Jesus

ADD tablemaker/tablemaker2.1.1.Linux_x86_64.tar.gz /

# Install tablemaker
RUN apt-get update
RUN cp -r tablemaker-2.1.1.Linux_x86_64/ /bin/
RUN rm -rf tablemaker-2.1.1.Linux_x86_64/

ENV PATH $PATH:/bin/tablemaker-2.1.1.Linux_x86_64
