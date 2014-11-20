NGSeasy (beta)
===================
[Funded by Biomedical Research Centre](http://core.brc.iop.kcl.ac.uk): http://core.brc.iop.kcl.ac.uk

Publication: pending

Authors: Stephen J Newhouse, Amos Folarin , Maximilian Kerz  
Release Version: 0.9.0  

****************
**A [Dockerized](https://www.docker.com/) and [Virtulaized](https://www.virtualbox.org/) ngs pipeline and tool-box.** 

**With NGSeasy you can now have full suite of NGS tools up and running on any high end workstation in an afternoon**

**Note:** NGSeasy is under **heavy development** and the code and docs evolve quickly.  

- **NGSeasy-v1.0 Full Production release will be available Dec 2014**    
- **NGSeasy updates every 12 months**
- **GUI in development**

****************

>NGSeasy is completely open source and we encourage interested folks to jump in and get involved in the dev with us.

****************

NGSeasy (Easy Analysis of Next Generation Sequencing)
=======================================================
We present **NGSeasy (Easy Analysis of Next Generation Sequencing)**, a flexible and easy-to-use NGS pipeline for automated alignment, quality control, variant calling and annotation. The pipeline allows users with minimal computational/bioinformatic skills to set up and run an NGS analysis on their own samples, in less than an afternoon, on any operating system (Windows, iOS or Linux) or infrastructure (workstation, cluster or cloud).

NGS pipelines typically utilize a large and varied range of software components and incur a substantial configuration burden during deployment which limits their portability to different computational environments. NGSeasy simplifies this by providing the pipeline components encapsulated in Docker™ containers and bundles in a wide choice of tools for each module. Each module of the pipeline represents one functional grouping of tools (e.g. sequence alignment, variant calling etc.).

Deploying the pipeline is as simple as pulling the container images from the public repository into any host running Docker. NGSeasy can be deployed on any medium to high-end workstation, high performance computer cluster and compute clouds (public/private cloud computing) - enabling instant access to elastic scalability without investment overheads for additional compute hardware and makes open and reproducible research straight forward for the greater scientific community.

### Advantages ###
- Easy to use for non-informaticians.  
- All run from a single config file that can be made in Excel.  
- User can select from mutiple aligners, variant callers and variant annotators
- No scary python, .yaml or .json files...just one simple Excel workbook saved as a textfile.  
- Just follow our simple set of instructions and NGS away!  
- Choice of aligners and variant callers and anntators  
- Allows reproducible research  
- Version controlled for auditing  
- Customisable  
- Easy to add new tools  
- If it's broke...we will fix it..
- Enforced naming convention and directory structures  
- Allows users to run "Bake Offs" between tools with ease  

We have adapted the current best practices from the Genome Analysis Toolkit (GATK, http://www.broadinstitute.org/gatk/guide/best-practices)  for processing raw alignments in SAM/BAM format and variant calling. The current workflow, has been optimised for Illumina platforms, but can easily be adapted for other sequencing platforms, with minimal effort.  

As the containers themselves can be run as executables with pre-specified cpu and RAM resources, the orchestration of the pipeline can be placed under the control of conventional load balancers if this mode is required.  

****
### Author Contact Details

Please contact us for help/guidance on using the beta release. 

- Amos Folarin <amosfolarin@gmail.com>  [@amosfolarin](https://twitter.com/amosfolarin?lang=en)   
<a href="http://www.linkedin.com/pub/amos-folarin/34/b06/978">
<img src="http://www.linkedin.com/img/webpromo/btn_viewmy_160x33.png" width="160" height="33" alt="View Amos's profile on LinkedIn">
</a>

- Stephen J Newhouse <stephen.j.newhouse@gmail.com> [@s_j_newhouse](https://twitter.com/s_j_newhouse?lang=en)  
<a href="http://uk.linkedin.com/pub/dr-stephen-newhouse/29/89a/11a">
<img src="http://www.linkedin.com/img/webpromo/btn_viewmy_160x33.png" width="160" height="33" alt="View Steve's profile on LinkedIn">
</a>

**Lets us know if you want other tools added to NGSeasy**

*Institution: NIHR Maudsley Biomedical Research Centre For Mental Health and Dementia Unit (Denmark Hill), at The Institute of Psychiatry, Psychology & Neuroscience (IoPPN), Kings College London* 

****************

Overview of the NGSeasy Pipeline Components
=============================================
The basic pipeline contains all the basic tools needed for manipulation and 
quality control of raw fastq files (ILLUMINA focused), SAM/BAM manipulation,
alignment, cleaning (based on GATK best practises [http://www.broadinstitute.org/gatk/guide/best-practices]) and first pass
variant discovery. Separate containers are provided for indepth variant annotation,
structural variant calling, basic reporting and visualisations.  

![ngsEASY](https://github.com/KHP-Informatics/ngs/blob/dev2/figs/ngsEASY_component_visualisation.png "Dockerized NGS Pipeline")


****************

The Full NGSeasy pipeline
=============================

The NGSeasy pipelines implement the following :-   

- **Quality control of raw fastq** files using **[FASTQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)**  

- **Read trimming** using **[TRIMMOMATIC](http://www.usadellab.org/cms/?page=trimmomatic)**.   

- **Alignment** using one of 
    - **[BWA](http://bio-bwa.sourceforge.net/)**  
    - **[STAMPY](http://www.well.ox.ac.uk/project-stampy)**   
    - **[NOVOALIGN](http://www.novocraft.com)**  
    - **[BOWTIE2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)**  
    - *[SNAP](http://snap.cs.berkeley.edu/): COMING SOON!*
    
- **SAM/BAM sorting and indexing** with **[SAMTOOLS](https://github.com/samtools/samtools)**.  

- **Read Group information added** using **[PICARDTOOLS](http://broadinstitute.github.io/picard/):[AddOrReplaceReadGroups](http://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups)** 

- **Duplicate marking** with **[PICARDTOOLS](http://broadinstitute.github.io/picard/):[MarkDuplicates](http://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates)**.  

>For academic users and/or commercial/clinical groups whom have paid for GATK licensing, the next steps are to perform   

- **Indel indel realignment and base quality score recalibration** using **[GATK](https://www.broadinstitute.org/gatk/)** built in tools :
    - **[GATK:RealignerTargetCreator](https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_indels_RealignerTargetCreator.php)** 
    - **[GATK:IndelRealigner](https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_indels_IndelRealigner.php)** 
    - **[GATK:BaseRecalibrator](https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php)** 

> For the non-GATK version    

- **Base quality score recalibration** using [BamUtil](http://genome.sph.umich.edu/wiki/BamUtil)  
    - **[BamUtil:recab](http://genome.sph.umich.edu/wiki/BamUtil:_recab)** 

- **Post alignment quality control and reporting** is performed usng a number of tools and custom scripts: 
    - **[SAMTOOLS:flagstats](https://github.com/samtools/samtools)**
    - **[BEDTOOLS:genomecov](https://github.com/arq5x/bedtools2)**
    - **[BEDTOOLS:bamtobed](https://github.com/arq5x/bedtools2)**
    - **[PICARDTOOLS:CollectMultipleMetrics](http://broadinstitute.github.io/picard/command-line-overview.html#CollectMultipleMetrics)**    
    - **[PICARDTOOLS:CollectAlignmentSummaryMetrics](http://broadinstitute.github.io/picard/command-line-overview.html#CollectAlignmentSummaryMetrics)**    
    - **[PICARDTOOLS:CollectWgsMetrics](http://broadinstitute.github.io/picard/command-line-overview.html#CollectWgsMetrics)**    
    - **[PICARDTOOLS:CollectTargetedPcrMetrics](http://broadinstitute.github.io/picard/command-line-overview.html#CollectTargetedPcrMetrics)** (coming soon)    

- **SNP and small INDEL** calling using one of the following or a combibation of these tools, if the `ensemble` method is called using **[bcbio.variation variant-ensemble](https://github.com/chapmanb/bcbio.variation)**
    - **[FREEBAYES](https://github.com/ekg/freebayes)** 
    - **[PLATYPUS](http://www.well.ox.ac.uk/platypus)** 
    - **[GATK:UnifiedGenotyper](https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_genotyper_UnifiedGenotyper.php)** 
    - **[GATK:HaplotypeCaller](https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php)** 

- **Structural Variant (CNV)** calling using one of the following or or a combibation of if the `ensemble` methods are called:- 
    - **[DELLY](https://github.com/tobiasrausch/delly)** 
    - **[LUMPY](https://github.com/arq5x/lumpy-sv/)**
    - **[cn.MOPS](http://www.bioinf.jku.at/software/cnmops/)**
    - **[m-HMM](https://www.stt.msu.edu/users/hengwang/mHMM.html)**
    - **[ExomeDepth](http://cran.r-project.org/web/packages/ExomeDepth/index.html)**
    - *[SLOPE](http://www-genepi.med.utah.edu/suppl/SLOPE/index.html) : COMING SOON!*

- **Variant annotation** using using one of the following or a combibation of if the `ensemble` methods are called. 
    - **[SnpEff](http://snpeff.sourceforge.net/)** 
    - **[ANNOVAR](http://www.openbioinformatics.org/annovar/)** 
    - **[VEP](http://www.ensembl.org/info/docs/tools/vep/index.html)**

- **Variant reporting** using custom scripts

**Note** Some of the later functions i.e. variant annotation and qc reporting are still in dev.  

****

**We highly recommed read trimming prior to alignment.**
We have noticed considerable speed-ups in alignmnet time and increased quality of SNP/INDEL calls using trimmed vs raw fastq.  

**Base quality score recalibration is also recommended.**  
As an alternative to GATK, we have added fucntionality for use of 
**[BamUtil](https://github.com/statgen/bamUtil):[recab](http://genome.sph.umich.edu/wiki/BamUtil:_recab)** 
for base quality score recalibration.

**Non-GATK users** 
- are encouraged to use aligners such as **[stampy](http://www.well.ox.ac.uk/project-stampy)** and **[novoalign](http://www.novocraft.com)** that perform base quality score recal on the fly.  
- are encouraged to use variant callers that perform local re-aligmnet around candidate sites to mitigate the need for the indel realignment stages.  
    - **[freebayes](https://github.com/ekg/freebayes)**
    - **[platypus](http://www.well.ox.ac.uk/platypus)**

********
### Coming Soon
- New Aligners:- [SNAP](http://snap.cs.berkeley.edu/), GSNAP, mr- and mrs-Fast,gem
- https://github.com/amplab/snap
- [SLOPE (CNV fo targetted NSG)] ((http://www.biomedcentral.com/1471-2164/12/184)) 
- Cancer Pipelines
- Annotation Pipelines and Databases
- Visualisation Pipelines
- Var Callers:- VarScan2
- SGE scripts and basic BASH scrips for running outside of Docker
- biobambam https://github.com/gt1/biobambam  
- bamaddrg https://github.com/ekg/bamaddrg  
- bamtools https://github.com/ekg/bamtools  


************

Dockerised NGSeasy
==========================
![docker](https://github.com/KHP-Informatics/ngs/blob/master/figs/Docker_container_engine_logo.png "Docker")  

The following section describes getting the Dockerised NGSeasy Pipeline(s) and Resources, project set up and running NGSeasy.

## 1. Install Docker

Follow the simple instructions in the links provided below  

- [Mac](https://docs.docker.com/installation/mac/)  
- [Windows](https://docs.docker.com/installation/windows/)
- [Ubuntu](https://docs.docker.com/installation/ubuntulinux/)

A full set of instructions for multiple operating systems are available on the [Docker website](https://docs.docker.com/installation/).

## 2. Get NGSeasy Automated build Container Images

All NGSeasy Docker images can be pulled down from **[compbio Docker Hub](https://hub.docker.com/u/compbio/)** or using the script [get_containers.sh](https://github.com/KHP-Informatics/ngs/blob/master/bin/get_containers.sh)

```bash
# get images
bash get_containers.sh v1.0
```

### Dockerised and Automated Builds ##

The following opensource tools are all provided as automated builds. 

| Tool | Build |
|-------------|----------------------|
|[ngseasy-base](https://registry.hub.docker.com/u/compbio/ngseasy-base/) | automated build |
|[fastqc](https://registry.hub.docker.com/u/compbio/ngseasy-fastqc) | automated build |
|[trimmomatic](https://registry.hub.docker.com/u/compbio/ngseasy-trimmomatic) | automated build |
|[bwa](https://registry.hub.docker.com/u/compbio/ngseasy-bwa) | automated build |
|[bowtie2](https://registry.hub.docker.com/u/compbio/ngseasybowtie-) | automated build |
|[picardtools](https://registry.hub.docker.com/u/compbio/ngseasy-picardtools) | automated build |
|[samtools](https://registry.hub.docker.com/u/compbio/ngseasy-samtools) | automated build |
|[freebayes](https://registry.hub.docker.com/u/compbio/ngseasy-freebayes/) | automated build |
|[bedtools](https://registry.hub.docker.com/u/compbio/ngseasy-bedtools/) | automated build |
|[bcbiovar](https://registry.hub.docker.com/u/compbio/ngseasy-bcbiovar/) | automated build |
|[delly](https://registry.hub.docker.com/u/compbio/ngseasy-delly) | automated build |
|[lumpy](https://registry.hub.docker.com/u/compbio/ngseasy-lumpy) | automated build |
|[cnmops](https://registry.hub.docker.com/u/compbio/ngseasy-cnmops) | automated build |
|[mhmm](https://registry.hub.docker.com/u/compbio/ngseasy-mhmm) | automated build |
|[exomedepth](https://registry.hub.docker.com/u/compbio/ngseasy-exomedepth) | automated build |
|[bamutil](https://registry.hub.docker.com/u/compbio/ngseasy-bamutil) | automated build |

samtools includes bcftools and htslib  

Its as easy as: - 
```bash
docker pull compbio/ngseasy-${TOOL}
```

******** 

## 3. Manually Build required NGSeasy Container Images

Currently we are not able to automatically build some of the tools in pre-built docker containers due to licensing restrictions. 

Some of the software has restrictions on use particularly for commercial 
purposes. Therefore if you wish to use this for commercial purposes, then you 
leagally have to approach the owners of the various components yourself!  

**Software composing the pipeline requiring registration:-**  

   * novoalign http://www.novocraft.com/  
   * Stampy http://www.well.ox.ac.uk/project-stampy  
   * Platypus http://www.well.ox.ac.uk/platypus  
   * GATK https://www.broadinstitute.org/gatk/  
   * ANNOVAR http://www.openbioinformatics.org/annovar/  

**These tools require manual download and registration with the proivder. For non-academics/commercial groups, you will need to pay for some of these tools.**

### Dockerised and Manual Builds ##

| Tool | Build |
|-------------|----------------------|
|[novoalign](https://github.com/KHP-Informatics/ngs/tree/master/containerized/ngs_docker_debian/ngseasy_novoalign) | manual build |
|[annovar](https://github.com/KHP-Informatics/ngs/tree/master/containerized/ngs_docker_debian/ngseasy_annovar) | manual build |
|[stampy](https://github.com/KHP-Informatics/ngs/tree/master/containerized/ngs_docker_debian/ngseasy_stampy) | manual build |
|[platypus](https://github.com/KHP-Informatics/ngs/tree/master/containerized/ngs_docker_debian/nsgeasy_platypus) | manual build |
|[gatk](https://github.com/KHP-Informatics/ngs/tree/master/containerized/ngs_docker_debian/ngseasy_gatk) | manual build |

Once you have paid/registered and downloaded the tool, we provide scripts and guidance for building these tools on your system.  

Its as easy as:-  
```{bash}
docker build -t compbio/ngseasy-${TOOL} .
```
******** 

## 4. Manually Build NGSeasy Variant Annotaion Container Images

The tools used for variant annotation use large databases and the docker images exceed 10GB. Therefore, the user should manually build these container images prior to running the NGS pipelines.
Docker build files ([Dockerfile](https://docs.docker.com/jsearch/?q=Dockerfile)) are available for 
- [Annovar](https://github.com/KHP-Informatics/ngs/tree/master/containerized/ngs_docker_debian/ngseasy_annovar/Dockerfile)  
- [VEP](https://github.com/KHP-Informatics/ngs/tree/master/containerized/ngs_docker_debian/ngseasy_vep/Dockerfile)   
- [snpEff](https://github.com/KHP-Informatics/ngs/tree/master/containerized/ngs_docker_debian/ngseasy_snpeff/Dockerfile)  

**Large Variant Annotation Container Images**

| Tool | Build |
|-------------|----------------------|
|[annovar](https://github.com/KHP-Informatics/ngs/tree/master/containerized/ngs_docker_debian/ngseasy_annovar) | manual build |
|[vep](https://github.com/KHP-Informatics/ngs/tree/master/containerized/ngs_docker_debian/ngseasy_vep) | manual build |
|[snpeff](https://github.com/KHP-Informatics/ngs/tree/master/containerized/ngs_docker_debian/ngseasy_snpeff) | manual build |


Its as easy as:-  
```{bash}
docker build -t compbio/ngseasy-${TOOL} .
```

**The user needs to manually build the following annotation tools:-   **

- Annovar
- VEP
- snpEff

**Note** Annovar requires user registration.

Once built on the user system, these container images can persist for as long as the user wants.

### Build VEP
```bash

cd /media/ngs_projects/nsgeasy/ngs/containerized/ngs_docker_debian/ngseasy_vep

sudo docker build -t compbio/ngseasy-vep:${VERSION} .
```

### Build Annovar

```bash
cd /media/ngs_projects/nsgeasy/ngs/containerized/ngs_docker_debian/ngseasy_annovar

sudo docker build -t compbio/ngseasy-annovar:${VERSION} .
```

### Build snpEff
```bash
cd /media/ngs_projects/nsgeasy/ngs/containerized/ngs_docker_debian/ngseasy_snpeff

sudo docker build -t compbio/ngseasy-snpeff:${VERSION} .
```

*******

## 5. Make NGS project directory
The user needs to make the relevent directory structures on their local machine before starting an NGS run. 

On our sysetm we typically set up a top-level driectory called **ngs_projects** within which we store output from all our individual NGS projects. 

Here we are working from local top level directory called **media/**, but this can really be any folder on your local system ie your home directory **~/${USER}**.  

Within this directory **media** we make the following folders: - 

```bash
ngs_projects  
|  
|__raw_fastq  
|__config_files  
|__reference_genomes_b37  
|__gatk_resources  
|__ngseasy
```

The following directories are obtained in **[6. Download NGSeasy Resources](https://github.com/KHP-Informatics/ngs/tree/master#6-download-ngseasy-resources)**
- **reference_genomes_b37**  
- **gatk_resources**  


```bash
# Move to media/
cd media

# make toplevel NGS folder
mkdir ngs_projects 

# fastq staging area
mkdir ngs_projects/fastq_raw 

# config files
mkdir ngs_projects/config_files 

# NGSeasy scripts
mkdir ngs_projects/ngseasy 
```

*****

## 6. Download NGSeasy Resources
Download the indexed reference genomes and example data for use with NGSeasy.

**NGSeasy Resources:-**  
- **reference_genomes_b37.tgz** b37 reference genomes indexed for use with all provided aligners (BWA, Bowtie2, Stampy, Novoalign) and annotation bed files for use with pipeline scripts
- **gatk_resources.tar.gz** gatk resources bundle
- **fastq_example.tgz** Example 75bp PE Illumina Whole Exome Sequence fastq data for **NA12878**
- Annotation Databases Coming in the next update 

**Download the data to the top level directory**

### FTP Details
- **ftp:**  159.92.120.21  
- **user:** compbio-public  
- **pwd:**  compbio-public  
- **port:** 21  

**Move to top level directory**
```bash
cd ngs_projects
```

**FTP NGSeasy Resources**
```bash
ftp 159.92.120.21
```

**mget NGSeasy Resources**
```bash
ftp> cd /Public/NGSeasy_Public_Resources
ftp> prompt off
ftp> mget *.gz
ftp> exit
```
I would recommend using a separate program like [FileZilla](https://filezilla-project.org/), which will make it much easier for you to set up and manage your file transfers

**Extract NGSeasy Resources**
```bash
# Extract resources
cd ngs_projects/
tar xvf gatk_resources.tar.gz; 
cd ngs_projects/gatk_resources
gunzip *

# Extract Reference Genomes
cd ngs_projects/
tar xvf reference_genomes_b37.tgz; 
cd ngs_projects/reference_genomes_b37
gunzip *
```

*****

## 7. Get NGSeasy Sripts
We then need to get the latest NGSeasy scripts from [GitHub](https://github.com/KHP-Informatics/ngs) . The user is required to download the scripts to the `ngseasy` directory

**move to the `ngseasy` directory**

```bash
cd /media/ngs_projects/nsgeasy
```

**clone the [ngs](https://github.com/KHP-Informatics/ngs) repository**

```bash
git clone https://github.com/KHP-Informatics/ngs.git
```

**add `nsgeasy/ngs/bin` to your system PATH**

```bash
export PATH=$PATH:/media/ngs_projects/nsgeasy/ngs/bin
```

**or add to global .bashrc**

```bash
echo "export PATH=$PATH:/media/ngs_projects/nsgeasy/ngs/bin" ~/.bashrc
source ~/.bashrc
```

**alternatively donwload the scripts from our [GitHub Release](https://github.com/KHP-Informatics/ngs)** 

### The NGSeasy Pipelines ##

| Pipeline             | Short Description    |
|----------------------|----------------------|
| [ngs_full_gatk](https://github.com/KHP-Informatics/ngs/blob/master/bin/ngs_full_gatk) | fastq to recalibrated bam to vcf using GATK  |
| [ngs_full_no_gatk](https://github.com/KHP-Informatics/ngs/blob/master/bin/ngs_full_no_gatk)    | fastq to recalibrated bam to vcf  |

gatk version includes indel realignment and base recalibration.  

Non-academics/commercial groups need to pay for GATK.  

Currently **ngs_full_gatk** pipeline is the most developed module.  

The **ngs_full_no_gatk** pipeline provides alternatives to processing with GATK. Here BamUtil:recab is used to recalibrate base quality scores and freebayes/platypus are the variant callers of choice.

****************

## 8. Set up NGS Project Working Directories

Running the script `ngseasy_initiate_project` ensures that all relevant directories are set up, and also enforces a clean structure to the NGS project.  

Within this we make a `raw_fastq` folder, where we temporarily store all the raw fastq files for each project. 
This folder acts as an initial stagging area for the raw fastq files. During the project set up, we copy/move project/sample related fastq files to their own specific directories.
Fastq files must have suffix and be gzipped: **_1.fq.gz** or **_2.fq.gz**  
furture version will allow any format  

Running `ngseasy_initiate_project` with the relevent configuration file, will set up the following directory structure for every project and sample within a project:-  

## NGS Project Directory 
```bash
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
	|__sample_id_1  
	|	|  
	|	|__fastq  
	|	|__tmp  
	|	|__alignments  
	|	|__vcf  
	|	|__reports  
	|	|__config_files  
	|
	|
	|__sample_id_n  
	|	|  
	|	|__fastq  
	|	|__tmp  
	|	|__alignments  
	|	|__vcf  
	|	|__reports  
	|	|__config_files  
```

**Running `ngseasy_initiate_project`**

```bash
ngseasy_initiate_project -c config.file.tsv -d /media/ngs_projects
```

****************

## 9. NGSeasy Project configuration file

In Excel make config file and save as [TAB] Delimited file with ``.tsv`` extenstion.  
See Example provided and [GoogleDoc](https://docs.google.com/spreadsheets/d/1kp1Nyw0x3zXqO2Wm2Z25ErJ0Z-Uoab8tjRPq9h4sonk/edit?usp=sharing). Remove the header from this file before running the pipeline. This sets up Information related to: Project Name, Sample Name, Library Type, Pipeline to call, NCPU.

The [config.file.tsv] should contain the following 15 columns for each sample to be run through a pipeline:- 

|Variable|type|Description|Options/Examples|
|--------|--------|--------|--------|
POJECT_ID|string|Project ID|Cancer|
SAMPLE_ID|string|Sample ID| T100|
FASTQ1|string|Raw fastq file name read 1| foo_1_fq.gz|
FASTQ2|string|Raw fastq file name read 1| foo_2_fq.gz|
PROJECT_DIR|string|Project Directory| /medida/ngs_projects |
DNA_PREP_LIBRARY_ID|string|DNA Libray Prep ID| Custom_Cancer |
NGS_PLATFORM|string|Platform Name| ILLUMINA |
NGS_TYPE|string|Experiment type| WGS/WEX/TGS/ |
BED_ANNO|string|Annotation Bed File|exons_b37.bed|
PIPELINE|string|NGSeasy Pipeline Script|ngs_full_gatk/ngs_full_no_gatk|
ALIGNER|string|Aligner| bwa/bowtie/stampy/novoalign|
VARCALLER|string|Variant Caller|ensemble/freebayes/platypus/UnifiedGenotyper/HaplotypeCaller|
GTMODEGATK|string|GATK Variant Caller Mode|EMIT_ALL_CONFIDENT_SITES/EMIT_VARIANTS_ONLY|
CLEANUP|string|Clean Up Files (TRUE/FALSE)|TRUE/FALSE|
NCPU|number|Number of cores to call|1..n|
VERSION|number|NGSeasy Version |v0.9/v1.0|


_coming soon_ options to add user email, specify non-gatk runs  

*************

## 10. Copy Project Fastq files to relevent Project/Sample Directories

```bash
ngseasy_initiate_fastq -c config.file.tsv -d /media/ngs_projects
```

****

## 11. Start the NGSeasy Volume Contaier

In the Docker container the project directory is mounted in `/home/pipeman/ngs_projects`

```bash
ngseasy_volumes_container -d /media/ngs_projects
```

****

## 12. Running a NGS full pipeline : from raw fastq to vcf calls

See https://github.com/KHP-Informatics/ngs/tree/dev2/bin for dev functions (Still working on these). 
Each of these will call a separate container and run a part of the NGS pipeline. Each step is usually 
dependent on the previous step(s) - in that they require certain data/input/output in the correct format 
and with the correct nameing conventions enforced by our pipeline to exist, before executing.

A full pipeline is set out below :-  


****

### Build Annotation Container Images



****

### Build Licensed and Registered Container Images

The user needs to manually build the following annotation tools:-   

- GATK
- Stampy
- Platypus
- Novoalign



****

### Set up project directories
```bash
ngseasy_initiate_project -c config.file.tsv -d /media/ngs_projects
```

### Copy fastq to project/sample directories
```bash
ngseasy_initiate_fastq -c config.file.tsv -d /media/ngs_projects
```

*** 

### A pipeline is called using :-

```bash
    ngseasy -c config.file.tsv -d /media/nsg_projects
```
    
In the config file we as to call the pipeline [ngs_full_gatk]. The pipeline is defined in the config file as [ngs_full_gatk].
Here [ngs_full_gatk] is a wrapper/fucntion for calling an NGS pipeline.  

Each function is a bash wrapper for an image/container(s) e.g:-  

**ngs_full_gatk**  

```bash
#!/bin/bash -x

#usage printing func
usage()
{
cat << EOF
  This script calls the NGSeasy pipeline ngs_full_gatk

  ARGUMENTS:
  -h      Flag: Show this help message
  -c      NGSeasy project and run configureation file
  -d      NGSeasy project directory

  EXAMPLE USAGE:
    
    ngseasy -c config.file.tsv -d project_directory

EOF
}

#get options for command line args
  while  getopts "hc:d:" opt
  do

      case ${opt} in
	  h)
	  usage #print help
	  exit 0
	  ;;
	  
	  c)
	  config_tsv=${OPTARG}
	  ;;

	  d)
	  project_directory=${OPTARG}
	  ;; 
      esac
  done

#check config file exists.
if [ ! -e "${config_tsv}" ] 
then
	    echo "ERROR :  ${config_tsv} does not exist "
	    usage;
	    exit 1;
fi

#check exists.
  if [ ! -d "${project_directory}" ] 
  then
	  echo " ERROR : ${project_directory} does not exist "
	  usage;
	  exit 1;
  fi

##################  
# start pipeline #
##################

ngseasy_fastqc  -c ${config_tsv} -d ${project_directory}

ngseasy_trimmomatic -c ${config_tsv} -d ${project_directory}

ngseasy_alignment -c ${config_tsv} -d ${project_directory}

ngseasy_addreadgroup -c ${config_tsv} -d ${project_directory}

ngseasy_markduplicates -c ${config_tsv} -d ${project_directory}

ngseasy_indel_realn -c ${config_tsv} -d ${project_directory}

ngseasy_base_recal -c ${config_tsv} -d ${project_directory}

ngseasy_filter_recalbam -c ${config_tsv} -d ${project_directory}

ngseasy_alignment_qc -c ${config_tsv} -d ${project_directory}
 
ngseasy_variant_calling -c ${config_tsv} -d ${project_directory}

# coming soon...
# ngseasy_filter_bam -c ${config_tsv} -d ${project_directory}
# ngseasy_cnv_calling -c ${config_tsv} -d ${project_directory}
# ngseasy_variant_filtering -c ${config_tsv} -d ${project_directory}
# ngseasy_variant_annotation -c ${config_tsv} -d ${project_directory}
# ngseasy_report -c ${config_tsv} -d ${project_directory}
```

****

Output suffixes 
===================

### Alignment Output
*.raw.sam  (WEX ~ 8GB)
*.raw.bam  
*.raw.bai  
*.sort.bam  (WEX ~ 3GB)
*.sort.bai  

****
### Addreadgroup
*.addrg.bam  
*.addrg.bai  
*.addrg.bam.bai  

****
### Dupemark
*.dupemk.bam  
*.dupemk.bai
*.dupemk.bam.bai  

***
### Indel realign
*.realn.bam  
*.realn.bai
*.realn.bam.bai  


***
### Base recal
*.recal.bam  (WEX ~ 4.4G)  
*.recal.bai  
*.recal.bam.bai  
*.realn.bam.BaseRecalibrator.table  
*.recal.bam.BaseRecalibrator.table  
*.recal.bam.BaseRecalibrator.BQSR.csv  

***

## Thresholds for Variant calling etc

For Freebayes and Platypus tools:-  

- We set min coverage to 10  
- Min mappinng quality to 20  
- Min base quality to 20

For GATK HaplotypeCaller (and UnifiedGenotyper)

```-stand_call_conf 30 -stand_emit_conf 10 -dcov 250 -minPruning 10```

Note: ```minPruning 10``` was added as many runs of HaplotypeCaller failed when using non-bwa aligend and GATK best practices cleaned BAMs. This fix sorted all problems out, and you really dont want dodgy variant calls...do you? Same goes for thresholds hard coded for use with Freebayes and Platypus.  
These setting all work well in our hands. Feel  free to edit the scripts to suit your needs.


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

### The Annoying thing about GATK!
This will break your runs if multiple calls try and access the file when the first call deletes it!  
```
WARN  11:05:27,577 RMDTrackBuilder - Index file /home/pipeman/gatk_resources/Mills_and_1000G_gold_standard.indels.b37.vcf.idx is out of date (index older than input file), deleting and updating the index file 
INFO  11:05:31,699 RMDTrackBuilder - Writing Tribble index to disk for file /home/pipeman/gatk_resources/Mills_and_1000G_gold_standard.indels.b37.vcf.idx 
```


## CNV tools to think about
EXCAVATOR: detecting copy number variants from whole-exome sequencing data @ http://genomebiology.com/2013/14/10/R120

>We developed a novel software tool, EXCAVATOR, for the detection of copy number variants (CNVs) from whole-exome sequencing data. EXCAVATOR combines a three-step normalization procedure with a novel heterogeneous hidden Markov model algorithm and a calling method that classifies genomic regions into five copy number states. We validate EXCAVATOR on three datasets and compare the results with three other methods. These analyses show that EXCAVATOR outperforms the other methods and is therefore a valuable tool for the investigation of CNVs in largescale projects, as well as in clinical research and diagnostics. EXCAVATOR is freely available at http://sourceforge.net/projects/excavatortool/ webcite.

*****

## Useful Links 

- https://bcbio.wordpress.com/  
- https://basecallbio.wordpress.com/2013/04/23/base-quality-score-rebinning/  
- https://github.com/statgen/bamUtil  
- http://genome.sph.umich.edu/wiki/BamUtil:_recab  
- https://github.com/chapmanb/bcbio.variation  
- http://plagnol-lab.blogspot.co.uk/2013/11/faq-and-clarifications-for-exomedepth-r.html
