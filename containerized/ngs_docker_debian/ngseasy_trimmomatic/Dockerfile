# Base image
FROM compbio/ngseasy-ubuntu-base:v1.0

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root

# Use baseimage-docker's init system.
CMD ["bash"]

# Set correct environment variables.
ENV HOME /root

# Update
RUN apt-get update -y

# + Trimmomatic
    RUN wget -O /tmp/Trimmomatic-0.32.zip http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip \
        && unzip /tmp/Trimmomatic-0.32.zip -d /usr/local/pipeline/ \
        && chown -R pipeman:ngsgroup /usr/local/pipeline/Trimmomatic-0.32 \
        && sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/pipeline/Trimmomatic-0.32/trimmomatic-0.32.jar' /home/pipeman/.bashrc \
        && sed -i '$aPATH=${PATH}:/usr/local/pipeline/Trimmomatic-0.32' /home/pipeman/.bashrc \
        && echo "alias ngsTrimfq='/usr/local/pipeline/Trimmomatic-0.32'" >>  /home/pipeman/.bashrc
       
#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 766 /usr/local/pipeline/***
RUN chown -R pipeman:ngsgroup /usr/local/pipeline

# Cleanup the temp dir
RUN rm -rf /tmp/*

# open ports private only
EXPOSE 80

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

#CMD and ENTRYPOINT
COPY ngseasy_hello.sh /home/pipeman/ngseasy_hello.sh
ENTRYPOINT ["/bin/bash"]
CMD ["/home/pipeman/ngseasy_hello.sh"]
