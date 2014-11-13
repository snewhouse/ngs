# Base image
FROM compbio/ngseasy-base:latest

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root
ENV DEBIAN_FRONTEND noninteractive

# Remain current
RUN apt-get update && \
  apt-get upgrade -y

# snpEff
RUN cd /usr/local/pipeline && \ 
  wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip && \
  unzip snpEff_latest_core.zip && \
  chown -R pipeman:ngsgroup /usr/local/pipeline/snpEff && \
  chmod -R 766 /usr/local/pipeline/snpEff/* && \
  sed -i '$aPATH=${PATH}:/usr/local/pipeline/snpEff' /home/pipeman/.bashrc && \
  echo "alias ngsSNPeff='/usr/local/pipeline/snpEff'" >>  /home/pipeman/.bashrc && \
  sed -i '$aPATH=${PATH}:/usr/local/pipeline/snpEff' /root/.bashrc

# Use baseimage-docker's bash.
CMD ["/bin/bash"]

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

