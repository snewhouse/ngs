NGSeasy
======================================

Getting and Setting up NGSeasy
***********************************

## Get NGSeasy & Resources

- git  
- ftp  
- docker  

## Local Machine Set Up

```
source /media/D/NGSeasy/nsg/ngseasy_scripts/ngseasy_functions
```

- Windows
- Mac

***********************************

Running a full pipeline from fastq to raw SNP/INDEL calls
-----------------------------------------------------------------------

## Step 1. Initiate NGSeasy Project
initiate_project -c /media/D/docker_ngs/ngseasy/config.file.tsv -d /media/D/docker_ngs/ngseasy/

## Step 2. Start Volumes Container
volumes_container -d /media/D/docker_ngs/ngseasy/

## Step 3.1 NGS Pipeline : Manual Calling (with gatk cleaning)

## Step 3.2 NGS Pipeline : Run Script (no gatk cleaning)

## Step 3.4 NGS Pipeline : Run Script (with gatk cleaning)

***********************************

Running aspects of the pipeline
-----------------------------------------------------------------------

***********************************

Runing NGSeasy Production Suite (AGENTS + DOCKER)
---------------------------------------------------

***********************************

Runing NGSeasy Production Suite (AGENTS + DOCKER)
---------------------------------------------------------

***********************************

Runing NGSeasy Annotation
-------------------------------

***********************************

