#########################################################################
# -- Author: Amos Folarin                                               #
# -- Email: amosfolarin@gmail.com                                       #
# -- Author: Stephen J Newhouse                                         #
# -- Email: stephen.j.newhouse@gmail.com                                #
# -- Organisation: KCL/SLaM                                             #
#########################################################################


#------------------------------------------------------------------------
# A simple script to acquire all the container images from DockerHub for
#  a given pipeline version.
#
# Either run as sudo or if you have docker configured to run without sudo
# then run as the regular user
# 
# Each major pipeline release has been given a tag, this file will enable 
# a record of which containers constituted the pipeline at that version
#------------------------------------------------------------------------

cat << EOF
# A simple script to acquire NGSeasy container images from DockerHub.
# Each major pipeline release has been given a tag, this file will enable
# a record of which containers constituted the pipeline at that version
# (see version tag in GitHub  https://github.com/KHP-Informatics/ngs )
#
# USAGE:
# Either run as sudo or if you have docker configured to run without sudo
# then run as the regular user and it will download the correct set of 
# container images
EOF


sudo docker pull compbio/ngseasy-fastx-toolkit:v0.9
sudo docker pull compbio/ngseasy-mhmm:v0.9
sudo docker pull compbio/ngseasy-exomedepth:v0.9
sudo docker pull compbio/ngseasy-cnmops:v0.9
sudo docker pull compbio/ngseasy-gsnap:v0.9
sudo docker pull compbio/ngseasy-gem:v0.9
sudo docker pull compbio/ngseasy-annovar:v0.9
sudo docker pull compbio/ngseasy-lumpy:v0.9
sudo docker pull compbio/ngseasy-delly:v0.9
sudo docker pull compbio/ngseasy-varscan2:v0.9
sudo docker pull compbio/ngseasy-freebayes:v0.9
sudo docker pull compbio/ngseasy-platypus:v0.9
sudo docker pull compbio/ngseasy-novoalign:v0.9
sudo docker pull compbio/ngseasy-stampy:v0.9
sudo docker pull compbio/ngseasy-bowtie2:v0.9
sudo docker pull compbio/ngseasy-vcftools:v0.9
sudo docker pull compbio/ngseasy-bedtools:v0.9
sudo docker pull compbio/ngseasy-seqtk:v0.9
sudo docker pull compbio/ngseasy-samtools:v0.9
sudo docker pull compbio/ngseasy-bwa:v0.9
sudo docker pull compbio/ngseasy-picardtools:v0.9
sudo docker pull compbio/ngseasy-gatk:v0.9
sudo docker pull compbio/ngseasy-trimmomatic:v0.9
sudo docker pull compbio/ngseasy-fastqc:v0.9
sudo docker pull compbio/ngseasy-base:wheezy



