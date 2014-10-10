# Base image
FROM compbio/ngseasy-base:wheezy

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root

# Update
RUN apt-get update -y

# + bowtie
RUN wget -O /tmp/bowtie2-2.2.3-linux-x86_64.zip http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.3/bowtie2-2.2.3-linux-x86_64.zip \
    && unzip /tmp/bowtie2-2.2.3-linux-x86_64.zip -d /usr/local/pipeline/ \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/bowtie2-2.2.3 \
    && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/bowtie2-2.2.3:/usr/local/pipeline/bowtie2-2.2.3/scripts' /home/pipeman/.bashrc \
    && echo "alias ngsBowtie2='/usr/local/pipeline/bowtie2-2.2.3'" >>  /home/pipeman/.bashrc

#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 766 /usr/local/pipeline/***
RUN chown -R pipeman:ngsgroup /usr/local/pipeline

# Cleanup the temp dir
RUN rm -rf /tmp/*

# open ports private only
EXPOSE 80

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN apt-get autoclean && apt-get autoremove -y && rm -rf /var/lib/{apt,dpkg,cache,log}/
