# NGSeasy Base Image

# FROM 
FROM debian:wheezy

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root

# Use baseimage-docker's bash.
CMD ["bash"]

## Remain current
RUN apt-get update -qq \
&& apt-get dist-upgrade -y

# basic dependencies
RUN apt-get update \
  && apt-get install -y \
  automake \
  ant \
  bash \
  binutils \
  bioperl \
  build-essential \
  bzip2 \
  c++11 \
  cdbs \
  cmake \
  cron \
  curl \
  dkms \
  dpkg-dev \
  g++ \
  gcc \
  gfortran \
  git \
  git-core \
  libblas-dev l\
  ibatlas-dev \
  libbz2-dev \
  liblzma-dev \
  libpcre3-dev \
  libreadline-dev \
  zlib1g-dev \  
  make \
  mercurial \
  openjdk-7-jdk \
  openjdk-7-jre \
  perl \
  php5-curl \
  python \
  sudo \
  subversion \
  tabix \
  tree \
  unzip \
  vim \
  wget
  

#------------------------------------------------------------------------
# USER SETUP
#------------------------------------------------------------------------
# Create a pipeline user:pipeman and group:ngsgroup
RUN useradd -m -s /bin/bash pipeman && cd /home/pipeman && echo "#bash config file for user pipeman" >> /home/pipeman/.bashrc
RUN groupadd ngsgroup
RUN usermod -aG ngsgroup pipeman
RUN usermod -aG sudo pipeman

#make pipeline install dirs
RUN mkdir /usr/local/pipeline && chown pipeman:ngsgroup /usr/local/pipeline && chmod 775 /usr/local/pipeline

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*



