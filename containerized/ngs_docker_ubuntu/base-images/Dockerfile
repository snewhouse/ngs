# FROM 
FROM ubuntu:latest

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root

# Use baseimage-docker's bash.
CMD ["bash"]

# Set correct environment variables.
ENV HOME /root

# basic dependencies
RUN apt-get update \
  && apt-get install -f -y \
  automake \
  bash \
  binutils \
  bioperl \
  build-essential \
  bzip2 \
  c++11 \
  cfitsio-dev \
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
  graphviz-dev \
  html2text \
  htmldoc \
  libavahi-compat-libdnssd-dev \
  libbam-dev \
  libboost-date-time-dev \
  libboost-filesystem-dev \
  libboost-iostreams-dev \
  libboost-program-options-dev \
  libboost-system-dev \
  libbz2-dev  \
  libcrypto++9 \
  libcrypto++-dev \
  libcurl3 \
  libcurl3-gnutls \
  libcurl4-openssl-dev \
  libfftw3-dev \
  libftgl-dev \
  libglew1.5-dev \
  libgsl0-dev \
  libjpeg-dev \
  libkrb5-dev \
  libldap2-dev \
  libmcrypt4 \
  libmcrypt-dev \
  libmpc3 \
  libmpc-dev \
  libmysqlclient-dev \
  libncurses5-dev \
  libpcre3-dev \
  libpng12-0 \
  libpng12-dev \
  libqt4-dev \
  libroot-dev \
  libssl1.0.0 \
  libssl-dev \
  libssl-dev \
  libx11-dev \
  libxext-dev \
  libxft-dev \
  libxml2 \
  libxml2-dev \
  libxml2-dev \
  libxpm-dev \
  libxslt1-dev \
  libxslt-dev \
  libzip-dev \
  linux-libc-dev \
  make \
  mercurial \
  openjdk-7-jdk \
  openjdk-7-jre \
  perl \
  php5-curl \
  python-dev \
  subversion \
  sysvbanner \
  tabix \
  tree \
  unzip \
  vim \
  wget \
  xlibmesa-glu-dev \
  zlib1g-dev;  


#------------------------------------------------------------------------
# USER SETUP
#------------------------------------------------------------------------
# Create a pipeline user:pipeman and group:ngsgroup
RUN useradd -m -s /bin/bash pipeman && cd /home/pipeman && echo "#bash config file for user pipeman" >> /home/pipeman/.bashrc
RUN groupadd ngsgroup
RUN usermod -G ngsgroup pipeman
RUN mkdir /usr/local/pipeline && chown pipeman:ngsgroup /usr/local/pipeline && chmod 775 /usr/local/pipeline

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*



