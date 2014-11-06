#------------------------------------------------#
# FASTQC TOOL: Dockerfile
#------------------------------------------------#

# Base image
FROM compbio/base-image:0.5

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root

# Update
RUN apt-get update -y

#-------------------------NGS-TOOL---------------------------------------
# + FastQC
    RUN wget -O /tmp/fastqc_v0.11.2.zip http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.2.zip \
        && unzip /tmp/fastqc_v0.11.2.zip -d /usr/local/pipeline/ \
        && chown -R pipeman:ngsgroup /usr/local/pipeline/FastQC \
        && sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/pipeline/FastQC/jbzip2-0.9.jar:/usr/local/pipeline/FastQC/sam-1.103.jar' /home/pipeman/.bashrc \
        && sed -i '$aPATH=${PATH}:/usr/local/pipeline/FastQC' /home/pipeman/.bashrc \
        && echo "alias ngsFastqc='/usr/local/pipeline/FastQC'" >>  /home/pipeman/.bashrc    

#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 766 /usr/local/pipeline/***
RUN chown -R pipeman:ngsgroup /usr/local/pipeline

#---------------------------------------------------------------------
#Cleanup the temp dir
RUN rm -rf /tmp/*

#open ports private only
EXPOSE 80

# CMD
CMD ["/bin/bash"]

#Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN apt-get autoclean && apt-get autoremove -y && rm -rf /var/lib/{apt,dpkg,cache,log}/

