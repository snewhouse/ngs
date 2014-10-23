#------------------------------------------------#
# bcbio.variation: Dockerfile
#------------------------------------------------#

# Base image
FROM compbio/ngseasy-base:wheezy

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root

# Update
RUN apt-get update -y && apt-get upgrade -y


#-------------------------NGS-TOOL---------------------------------------
RUN mkdir /usr/local/pipeline/bcbio  && \
  cd /usr/local/pipeline/bcbio && \
  wget https://github.com/chapmanb/bcbio.variation/releases/download/v0.1.7/bcbio.variation-0.1.7-standalone.jar && \
  sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/pipeline/bcbio/bcbio.variation-0.1.7-standalone.jar' /home/pipeman/.bashrc && \
  sed -i '$aPATH=${PATH}:/usr/local/pipeline/bcbio/' /home/pipeman/.bashrc && \
  sed -i '$aPATH=${PATH}:/usr/local/pipeline/bcbio/' ~/.bashrc && \
  echo "alias ngsBCBIOVAR'/usr/local/pipeline/bcbio/'" >>  /home/pipeman/.bashrc

#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 777 /usr/local/pipeline
RUN chown -R pipeman:ngsgroup /usr/local/pipeline

#---------------------------------------------------------------------
#Cleanup the temp dir
RUN rm -rf /tmp/*

#open ports private only
EXPOSE 80

#Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN apt-get autoclean && apt-get autoremove -y && rm -rf /var/lib/{apt,dpkg,cache,log}/





