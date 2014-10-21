NGSeasy
==========================

....Under Construction...

NGSeasy Project Set up
--------------------------

## Step 0. Set up project directories

Fastq files must have suffix and be gzipped: **_1.fq.gz** or **_2.fq.gz**  
furture version will allow any format  

```{bash}
#--------------------------------#
# make top level dirs 
#--------------------------------#
cd media
mkdir ngs_projects
mkdir ngs_projects/fastq_raw
mkdir ngs_projects/config_files
```
****

## Step 1. Set up project configuration file

In Excel make config file and save as [TAB] Delimited file with ``.tsv`` extenstion.  
See Example provided and [GoogleDoc](https://docs.google.com/spreadsheets/d/1kp1Nyw0x3zXqO2Wm2Z25ErJ0Z-Uoab8tjRPq9h4sonk/edit?usp=sharing). Remove the header from this file before running the pipeline. This sets up Information related to: Project Name, Sample Name, Library Type, Pipeline to call, NCPU.

The [config.file.tsv] should contain the following 15 columns for each sample to be run through a pipeline:- 

|Variable|type|Description|
|--------|--------|--------|
POJECT_ID|string|Project ID|
SAMPLE_ID|string|Sample ID|
FASTQ1|string|Raw fastq file name read 1|
FASTQ2|string|Raw fastq file name read 1|
PROJECT_DIR|string|Project Directory|
DNA_PREP_LIBRARY_ID|string|DNA Libray Prep ID|
NGS_PLATFORM|string|Platform Name|
NGS_TYPE|string|Experiment type|
BED_ANNO|string|Annotation Bed File
PIPELINE|string|NGSeasy Pipeline Script|
ALIGNER|string|Aligner|
VARCALLER|string|Variant Caller|
GTMODEGATK|string|GATK Variant Caller Mode|
CLEANUP|string|Clean Up Files (TRUE/FALSE)|
NCPU|number|Number of cores to call|

****

## Step 2. Initiate the Project
The user needs to make the relevent directory structure on their local machine. Running this script ensures that all relevant directories are set up, ans also enforces a clean structure to the NGS project.  

On our sysetm we typically set up a top-level driectory called `ngs_projects` within which we store output from all our individual NGS projects. Within this we make a `raw_fastq` folder, where we temporarily store all the raw fastq files for each project. This folder acts as an initial stagging area for the raw fastq files. During the project set up, we copy/move project/sample related fastq files to their own specific directories.

Running `ngseasy_initiate_project` with the relevent configuration file, will set up the following directory structure for every project and sample within a project:-  

### NGS Project Directory 
```
.
ngs_projects  
|  
|__raw_fastq  
|__config_files  
|__reference_genomes_b37  
|__gatk_resources  
|__ngseasy
|
|__ project_id  
	|  
	|__run_logs  
	|__config_files  
	|__project_vcfs  
	|__project_bams  
	|__project_reports  
	|
	|__sample_id  
		|  
		|__fastq  
		|__tmp  
		|__alignments  
		|__vcf  
		|__reports  
		|__config_files  

```
### Running **ngseasy_initiate_project**

```{bash}
ngseasy_initiate_project -c config.file.tsv -d /media/ngs_projects
```
****

## Step 3. Copy Project Fastq files to relevent Project/Sample Directories

```{bash}
ngseasy_initiate_fastq -c config.file.tsv -d /media/ngs_projects
```

****

## Step 4. Start the NGSeasy Volume Contaier

In the Docker container the project directory is mounted in `/home/pipeman/ngs_projects`

```{bash}
ngseasy_volumes_container -d /media/ngs_projects
```

## Summary

```{bash}
#make top level dirs 

mkdir ngs_projects
mkdir ngs_projects/fastq_raw
mkdir ngs_projects/config_files

#copy/download raw fastq file to [ngs_projects/fastq_raw]

#set up project specific configuration file [config.file.tsv]

# Set up NGSeasy
ngseasy_initiate_project -c config.file.tsv -d /media/ngs_projects
ngseasy_initiate_fastq -c config.file.tsv -d /media/ngs_projects
ngseasy_volumes_container -d /media/ngs_projects
```

** You are now ready to run **

****

## Running a full pipeline : from raw fastq to vcf calls

See https://github.com/KHP-Informatics/ngs/tree/dev2/bin for dev functions (Still workin on these). Each of these will call a separate container and run a part of the NGS pipeline. Each step is usually dependent on the previous step(s) - in that they require certain data/input/output in the correct format and with the correct nameing conventions enforced by our pipeline to exist, before executing.

A full pipeline is set out below :-  

```{bash}

#--------------------------------#
# make top level dirs 
#--------------------------------#
cd media
mkdir ngs_projects
mkdir ngs_projects/fastq_raw # fastq staging area
mkdir ngs_projects/config_files # config files
mkdir ngs_projects/humandb # for annovar databses

#get NGSeasy resources
# sftp From ........copy data to and extract
cd ngs_projects
sftp....

tar xvf gatk_resources.tar.gz; gunzip *
tar xvf reference_genomes_b37.tgz; gunzip *


#--------------------------------#
# get and PATH nsgeasy scripts
#--------------------------------#

cd ngs_projects/nsgeasy
git clone https://github.com/KHP-Informatics/ngs.git
git checkout dev2
export PATH=$PATH:/media/ngs_projects/nsgeasy/ngs/bin
ln -s /media/ngs_projects/nsgeasy/ngs/ngeasy_dev/bin /media/ngs_projects/ngseasy_scripts

#to do [get_annovar_humandb]

#get images
bash get_containers.sh

#------------------------------------------------#
# to be run outside of docker and before ngseasy #
#------------------------------------------------#

ngseasy_initiate_project -c config.file.tsv -d /media/ngs_projects

ngseasy_initiate_fastq -c config.file.tsv -d /media/ngs_projects

ngseasy_get_annovar_db -d /media/ngs_projects/humandb

#--------------------#
# NGSEASY Dockerised #
#--------------------#

# A pipeline is called using :-

    ngseasy -c config.file.tsv -d /media/nsg_projects

# in the config file we as to call the pipeline [full]
# here [ngs_full_gatk] is a wrapper/fucntion for calling the pipeline

ngs_full_gatk() { 


#step through a full pipeline

#1 FastQC raw data
ngseasy_fastqc -c ${config_tsv} -d ${project_directory} -t 1;

#2 Trimm Fastq files
ngseasy_fastq_trimm -c ${config_tsv} -d ${project_directory};

#3 FastQC trimmed file
ngseasy_fastqqc -c ${config_tsv} -d ${project_directory} -t 0;

#4 Run Alignment
ngseasy_alignment -c ${config_tsv} -d ${project_directory};

#5 Add read groups ie sample info to BAM file
nsgeasy_add_readgroups -c ${config_tsv} -d ${project_directory};

#5 Mark Duplicates
ngseasy_mark_dupes -c ${config_tsv} -d ${project_directory};

#6 Indel Realignment
nsgeasy_indel_realignment -c ${config_tsv} -d ${project_directory};

#7 Base Recalibration
nsgeasy_base_recalibration -c ${config_tsv} -d ${project_directory};

#8 Alignment QC reports
ngseasy_alignment_qc -c ${config_tsv} -d ${project_directory};

#9 Call Vairants : SNPS and small Indels
nsgeasy_var_callers -c ${config_tsv} -d ${project_directory};

#10 CNV Calling
nsgeasy_cnv_callers -c ${config_tsv} -d ${project_directory};

#11 Variant Annotation: Annovar
nsgeasy_variant_annotation -c ${config_tsv} -d ${project_directory};

#TODO
#variant filter
#nsgeasy_variant_combine -c ${config_tsv} -d ${project_directory}
#nsgeasy_variant_report -c ${config_tsv} -d ${project_directory}


}

recommend full : trimmed aln gatk filtered and ensemble calls (multi SNP/INDELS/SV callers) base recalibration
if not novoalign then stampy (bwa with stampy)


```
****

To add :
- getting the pipeline and setting up resources data
- not all steps need config?
- pipeline option need to be set how? list of steps, specified full, full_no_gatk, var_call_only, cnv_call_only, qc_reports??

****
### Alignment Output
*.raw.sam
*.raw.bam
*.raw.bai
*.sort.bam
*.sort.bai

****
### Addreadgroup
*.addrg.bam
*.addrg.bam.bai

****
### Dupemark
*.dupemk.bam
*.dupemk.bam.bai


****
blah blah blah

****

## Gottchas

**bin/bash -c**

- need to add ```/bin/bash -c ${COMMAND}``` when software require ```>``` redirect to some output

example below for bwa:-  

```
  sudo docker run \
  -P \
  --name sam2bam_${SAMPLE_ID} \
  --volumes-from volumes_container \
  -t compbio/ngseasy-samtools:v0.9 /bin/bash -c \
  "/usr/local/pipeline/samtools/samtools view -bhS ${SOUTDocker}/alignments/${BAM_PREFIX}.raw.bwa.sam > ${SOUTDocker}/alignments/${BAM_PREFIX}.raw.bwa.bam"
  ```

runnig this without ```/bin/bash -c``` breaks. The ```>``` is called outside of the container