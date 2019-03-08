FROM python:2.7.15

MAINTAINER Filipe Jesus

ADD miso/complete_run.sh /

RUN apt-get update
RUN pip install misopy==0.5.4
RUN pip install bsddb3

RUN apt-get install bedtools

#Download samtools from GitHub
WORKDIR /bin
RUN wget http://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2
RUN tar --bzip2 -xf samtools-1.2.tar.bz2
WORKDIR /bin/samtools-1.2
RUN make
RUN rm /bin/samtools-1.2.tar.bz2

ENV PATH $PATH:/bin/samtools-1.2

RUN cd /usr/local/lib/python2.7/site-packages/misopy/ && mkdir settings && cd settings && touch /miso_settings.txt && echo "[data]" >> miso_settings.txt && echo "filter_results = False" >> miso_settings.txt && echo "min_event_reads = 20" >> miso_settings.txt && echo "strand = fr-unstranded" >> miso_settings.txt && echo "" >> miso_settings.txt && echo "[cluster]" >> miso_settings.txt && echo "cluster_command = qsub" >> miso_settings.txt && echo "" >> miso_settings.txt && echo "[sampler]" >> miso_settings.txt && echo "burn_in = 500" >> miso_settings.txt && echo "lag = 10" >> miso_settings.txt && echo "num_iters = 5000" >> miso_settings.txt && echo "num_processors = 4" >> miso_settings.txt && cd /
RUN chmod a=wrx /usr/local/lib/python2.7/site-packages/misopy/settings/miso_settings.txt

WORKDIR /

RUN apt-get clean
