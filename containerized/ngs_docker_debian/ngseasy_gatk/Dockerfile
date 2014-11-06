#------------------------------------------------#
# GATK: Dockerfile
#------------------------------------------------#

# Base image
FROM compbio/base-image

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root

# Update
RUN apt-get update -y

#set java 
ENV JAVA_HOME /usr/lib/jvm/java-1.7.0-openjdk-amd64
ENV PATH:${PATH}:${JAVA_HOME}

#-------------------------NGS-TOOL---------------------------------------
ADD GenomeAnalysisTK-3.2-2.tar.bz2 /usr/local/pipeline/GenomeAnalysisTK-3.2-2

RUN sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/pipeline/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar' /home/pipeman/.bashrc && \
    sed -i '$aPATH=${PATH}:/usr/local/pipeline/GenomeAnalysisTK-3.2-2' /home/pipeman/.bashrc && \
    echo "alias ngsGATK='/usr/local/pipeline/GenomeAnalysisTK-3.2-2'" >>  /home/pipeman/.bashrc && \

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
