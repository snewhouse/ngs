NGSeasy
==========================

....Under Construction...

NGSeasy Project Set up
--------------------------

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


## Step 2. Initiate the Project
The user needs to make the relevent directory structure on their local machine. Running this script ensures that all relevant directories are set up, ans also enforces a clean structure to the NGS project.  

On our sysetm we typically set up a top-level driectory called `ngs_projects` within which we store output from all our individual NGS projects. Within this we make a `raw_fastq` folder, where we temporarily store all the raw fastq files for each project. This folder acts as an initial stagging area for the raw fastq files. During the project set up, we copy/move project/sample related fastq files to their own specific directories.

Running `ngseasy_initiate_project` with the relevent configuration file, will set up the following directory structure for every project and sample within a project:-  

```
.
ngs_projects  
|  
|__raw_fastq
|
|__ project_id  
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
Running **ngseasy_initiate_project**

```{bash}
ngseasy_initiate_project -c config.file.tsv -d /media/ngs_projects
```

## Step 3. Copy Project Fastq files to relevent Project/Sample Directories

```{bash}
ngseasy_initiate_fastq -c config.file.tsv -d /media/ngs_projects
```

## Step 4. Start the NGSeasy Volume Contaier

In the Docker container the project directory is mounted in `/home/pipeman/ngs_projects`

```{bash}
ngseasy_volumes_container -d /media/ngs_projects
```

## Summary

```{bash}
#make top level dirs 
mkdir ngs_projects
mkdir raw_fastq

#copy/download raw fastq file to [ngs_projects/raw_fastq]

#set up project specific configuration file [config.file.tsv]

# Set up NGSeasy
ngseasy_initiate_project -c config.file.tsv -d /media/ngs_projects
ngseasy_initiate_fastq -c config.file.tsv -d /media/ngs_projects
ngseasy_volumes_container -d /media/ngs_projects
```

** You are now ready to run **

****

## Running a full pipeline : from raw fastq to vcf calls

See https://github.com/KHP-Informatics/ngs/tree/dev2/ngseasy_dev for dev functions (Still workin on these). Each of these will call a separate container and run a part of the NGS pipeline. A full pipeline is set out below :-  

```{bash}

ngseasy_initiate_project -c config.file.tsv -d /media/ngs_projects

ngseasy_initiate_fastq -c config.file.tsv -d /media/ngs_projects

ngseasy_volumes_container -d /media/ngs_projects

ngseasy_fastq_qc -c config.file.tsv -d /media/ngs_projects

ngseasy_fastq_trimm -c config.file.tsv -d /media/ngs_projects

ngseasy_alignment -c config.file.tsv -d /media/ngs_projects

ngseasy_mark_dupes -c config.file.tsv -d /media/ngs_projects

nsgeasy_add_readgroups -c config.file.tsv -d /media/ngs_projects

nsgeasy_indel_realignment -c config.file.tsv -d /media/ngs_projects

nsgeasy_base_recalibration -c config.file.tsv -d /media/ngs_projects

ngseasy_alignment_qc -c config.file.tsv -d /media/ngs_projects

nsgeasy_var_callers -c config.file.tsv -d /media/ngs_projects

nsgeasy_cnv_callers -c config.file.tsv -d /media/ngs_projects

nsgeasy_variant_annotation -c config.file.tsv -d /media/ngs_projects

nsgeasy_variant_combine -c config.file.tsv -d /media/ngs_projects

nsgeasy_variant_report -c config.file.tsv -d /media/ngs_projects

```
****

To add :
- getting the pipeline and setting up resources data
- not all steps need config?
- pipeline option need to be set how? list of steps, specified full, full_no_gatk, var_call_only, cnv_call_only, qc_reports??