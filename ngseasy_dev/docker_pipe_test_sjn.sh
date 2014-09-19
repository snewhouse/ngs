#!/bin/bash
# Collection of NGSeasy Functions
# Stephen Newhouse <stephen.j.newhouse@gmail.com>
# Version 0.9.0

##--------------------------------------------------##
## create project directories
##--------------------------------------------------##

initiate_project () {

# NSGeasy Project Initiation
# Makes Project and Sample Directories needed for NGSesy
# Stephen Newhouse <stephen.j.newhouse@gmail.com>
# Amos Folarin <afolarin@gmail.com>

#usage printing func
usage()
{
cat << EOF
This script sets up the NGSeasy Project directories
See NGSEasy containerized instructions.

ARGUMENTS:
-h      Flag: Show this help message
-c      Config pipeline file
-d      Base directory for (fastq_raw, reference_genomes_b37, gatk_resources, ngs_projects, ngseasy_scripts)
EXAMPLE USAGE:
initiate_project -c config.file.tsv -d /media/D/docker_ngs/ngseasy/
EOF
}


#get options for command line args
while  getopts "h:c:d:" opt
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
        host_vol_dir=${OPTARG}
        ;;

    esac
done

#check exists.
  if [[ ! -e ${host_vol_dir} ]] 
  then
	  echo " ${host_vol_dir} does not exist "
	  usage;
	  exit 1;
  fi
#check exists.
  if [[ ! -e ${config_tsv} ]] 
  then
	  echo " ${config_tsv} does not exist "
	  usage;
	  exit 1;
  fi
  

## Print to screen input 
echo  " NGS config_file  : [${config_file}]"
echo  " NGS host_vol_dir : [${host_vol_di}]"

echo  " Checking [${config_file}] format"

  numfeilds=`awk '{print NF}' ${config_file} | sort -gr | sed '1q'`
  numsamples=`wc -l ${config_file} | awk '{print $1}'`
echo  " Number of samples : [$numsamples]"
    
  if [ "$numfeilds" -ne "15" ]
    then
    echo  " WARNING! Number of fields not equal to 15 for one of your samples You appear to have missing data in [${config_file}]"
    echo  " Check  [${config_file}] format"
    echo  " The [config.file] should contain the following 15 columns for each sample:-"
    echo  " 
      POJECT_ID
      SAMPLE_ID
      FASTQ1
      FASTQ2
      PROJECT_DIR
      DNA_PREP_LIBRARY_ID
      NGS_PLATFORM
      NGS_TYPE
      BED_ANNO
      PIPELINE
      ALIGNER
      VARCALLER
      GTMODEGATK
      CLEANUP
      NCPU"
    echo  " Data are [TAB] delimited and each line should end in a [NEW LINE]"
      exit 1
    else
   echo  " All Good...Proceeding to read : [${config_file}]"
  fi

# begin reading config file line by line
  
while read -r f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15
do

# set varibales  
  POJECT_ID=$f1
  SAMPLE_ID=$f2
  FASTQ1=$f3
  FASTQ2=$f4
  PROJECT_DIR=$f5
  DNA_PREP_LIBRARY_ID=$f6
  NGS_PLATFORM=$f7
  NGS_TYPE=$f8
  BED_ANNO=$f9
  PIPELINE=$f10
  ALIGNER=$f11
  VARCALLER=$f12
  GTMODEGATK=$f13
  CLEANUP=$f14
  NCPU=$f15
  
echo  "#------------------------------------------------------------#"
echo  " Checking Project/Sample Directory Structure"
echo  "#------------------------------------------------------------#"
     echo " Checking Project Directory "   
  if [ ! -d $PROJECT_DIR/${POJECT_ID} ]
    then
    echo  " Project Directory Does not exist" 
    echo  " Setting up correct Project Directory"
    echo ""
    echo  " Making Directory : [$PROJECT_DIR/${POJECT_ID}/]"
    echo  " Making Directory : [$PROJECT_DIR/${POJECT_ID}/config_files]"
    echo  " Making Directory : [$PROJECT_DIR/${POJECT_ID}/cohort_vcfs]"
    
        mkdir ${PROJECT_DIR}/${POJECT_ID}/
        mkdir ${PROJECT_DIR}/${POJECT_ID}/config_files
        mkdir ${PROJECT_DIR}/${POJECT_ID}/cohort_vcfs
        
  else
    echo  " Project Directory Exists"
  fi

     echo " Checking Sample Directory " 
  if [ ! -d $PROJECT_DIR/${POJECT_ID}/$SAMPLE_ID ]  
    then
    echo " Sample Directory Does not exist"
    echo " Setting up correct Sample Directory"
    echo ""
        mkdir ${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}
        mkdir ${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/fastq
        mkdir ${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/tmp
        mkdir ${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/alignments
        mkdir ${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/vcf
        mkdir ${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/reports
        mkdir ${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/config_files
        
    echo  " Making Sample Directory : [${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}]"
    echo  " Making Sample Directory : [${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/fastq]"
    echo  " Making Sample Directory : [${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/tmp]"
    echo  " Making Sample Directory : [${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/alignments]"
    echo  " Making Sample Directory : [${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/vcf]"
    echo  " Making Sample Directory : [${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/reports]"
    echo  " Making Sample Directory : [${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/config_files]"
  
  else
    echo  " Sample Directory Exists "
  fi
  
echo ""
echo  "#------------------------------------------------------------#"
echo  " NGSeasy Project/Sample Settings "
echo  "#------------------------------------------------------------#"
DATE=`date +"%d%m%y"`
echo  " Run Date : [$DATE]"
echo  " Project_Id : [$POJECT_ID]" 
echo  " Project Directory : [$PROJECT_DIR]"
echo  " Platform : [$NGS_PLATFORM]"
echo  " NGS Type : [$NGS_TYPE]"
echo  " Clean up TRUE/FALE : [$CLEANUP]"
echo  " Number of cpu : [$NCPU]"
echo  " Output Directory : [${PROJECT_DIR}/${POJECT_ID}]"
echo  " Saving Project configuration settings to : [${PROJECT_DIR}/${POJECT_ID}/config_files/project.config.file.${DATE}]"
echo  "$DATE $f1 $f2 $f3 $f4 $f5 $f6 $f7 $f8 $f9 $f10 $f11 $f12 $f13 $f14 $f15 ${SAMPLE_ID}.${NGS_TYPE}.${NGS_PLATFORM}.${ALIGNER} ${PROJECT_DIR}/${POJECT_ID}" >> ${PROJECT_DIR}/${POJECT_ID}/config_files/project.config.file.${DATE};
  
done < ${config_file}
}

##--------------------------------------------------##
## create and run data_volumes container 
##--------------------------------------------------##

volumes_container() {

#usage printing func
  usage()
  {
  cat << EOF
  This script sets up the docker volumes container:
  See NGSEasy containerized instructions.

  ARGUMENTS:
  -h      Flag: Show this help message
  -d      Base directory for (fastq_raw, reference_genomes_b37, gatk_resources, ngs_projects, ngseasy_scripts)
  EXAMPLE USAGE:
  volumes_container -d /media/D/docker_ngs/ngseasy/
EOF 
}

#get options for command line args
  while  getopts "h:d:" opt
  do

      case ${opt} in
	  h)
	  usage #print help
	  exit 0
	  ;;
	  
	  d)
	  host_vol_dir=${OPTARG}
	  ;;

      esac
  done

#check exists.
  if [[ ! -e ${host_vol_dir} ]] 
  then
	  echo " ${host_vol_dir} does not exist "
	  usage;
	  exit 1;
  fi

#run docker image  
  sudo docker run \
  -d \
  -P \
  -v ${host_vol_dir}/fastq_raw:/home/pipeman/fastq_raw \
  -v ${host_vol_dir}/reference_genomes_b37:/home/pipeman/reference_genomes_b37 \
  -v ${host_vol_dir}/gatk_resources:/home/pipeman/gatk_resources \
  -v ${host_vol_dir}/ngs_projects:/home/pipeman/ngs_projects \
  -v ${host_vol_dir}/ngseasy_scripts:/home/pipeman/ngseasy_scripts \
  --name data_volumes \
  -t compbio/ngseasy-base:v1.0
}

 

##--------------------------------------------------##
## NGSeasy: FastQC Pre-Alignment
##--------------------------------------------------##

fastqc_pre () {
#usage printing func
  usage()
  {
  cat << EOF
  This script sets up the NGSeasy docker fastqc container:
  See NGSEasy containerized instructions.

  ARGUMENTS:
  -h      Flag: Show this help message
  -c      NGSeasy project and run configureation file
  EXAMPLE USAGE:
  fastqc_pre -c config.file.tsv
EOF 
}

#get options for command line args
  while  getopts "h:c:" opt
  do

      case ${opt} in
	  h)
	  usage #print help
	  exit 0
	  ;;
	  
	  c)
	  config_tsv=${OPTARG}
	  ;;

      esac
  done

#check exists.
  if [[ ! -e ${config_tsv} ]] 
  then
	  echo " ${config_tsv} does not exist "
	  usage;
	  exit 1;
  fi

#run compbio/ngseasy-fastq
sudo docker run \
-d \
-P \
--name fastqc_pre \
--volumes-from volumes_container \
-t compbio/ngseasy-fastqc:v1.0 /sbin/my_init -- /bin/bash /home/pipeman/ngseasy_scripts/run_ngseasy_fastqc_pre_trimm.sh /home/pipeman/ngs_projects/${config_tsv}

}

##--------------------------------------------------##
## NGSeasy: FastQC Trimmomatic
##--------------------------------------------------##

fastq_trimm () {
#usage printing func
  usage()
  {
  cat << EOF
  This script sets up the NGSeasy docker fastqc Trimmomatic container:
  See NGSEasy containerized instructions.

  ARGUMENTS:
  -h      Flag: Show this help message
  -c      NGSeasy project and run configureation file
  EXAMPLE USAGE:
  fastq_trimm -c config.file.tsv
EOF 
}

#get options for command line args
  while  getopts "h:c:" opt
  do

      case ${opt} in
	  h)
	  usage #print help
	  exit 0
	  ;;
	  
	  c)
	  config_tsv=${OPTARG}
	  ;;

      esac
  done

#check exists.
  if [[ ! -e ${config_tsv} ]] 
  then
	  echo " ${config_tsv} does not exist "
	  usage;
	  exit 1;
  fi

sudo docker run \
-d \
-P \
--name fastq_trimm \
--volumes-from volumes_container \
-t compbio/ngseasy-fastqc:v1.0 /sbin/my_init -- /bin/bash /home/pipeman/ngseasy_scripts/run_ngseasy_trimmomatic.sh /home/pipeman/ngs_projects/${config_tsv}

}


##--------------------------------------------------##
## NGSeasy: FastQC Post-Trimmomatic
##--------------------------------------------------##

fastqc_post () {
#usage printing func
  usage()
  {
  cat << EOF
  This script sets up the NGSeasy docker fastqc container:
  See NGSEasy containerized instructions.

  ARGUMENTS:
  -h      Flag: Show this help message
  -c      NGSeasy project and run configureation file
  EXAMPLE USAGE:
  fastqc_post -c config.file.tsv
EOF 
}

#get options for command line args
  while  getopts "h:c:" opt
  do

      case ${opt} in
	  h)
	  usage #print help
	  exit 0
	  ;;
	  
	  c)
	  config_tsv=${OPTARG}
	  ;;

      esac
  done

#check exists.
  if [[ ! -e ${config_tsv} ]] 
  then
	  echo " ${config_tsv} does not exist "
	  usage;
	  exit 1;
  fi

#run compbio/ngseasy-fastq
sudo docker run \
-d \
-P \
--name fastqc_post \
--volumes-from volumes_container \
-t compbio/ngseasy-fastqc:v1.0 /sbin/my_init -- /bin/bash /home/pipeman/ngseasy_scripts/run_ngseasy_fastqc_post_trimm.sh /home/pipeman/ngs_projects/${config_tsv}

}

##--------------------------------------------------##
## NGSeasy: aln_novo
##--------------------------------------------------##

aln_novo () {
#usage printing func
  usage()
  {
  cat << EOF
  This script sets up the NGSeasy docker novoalign container:
  See NGSEasy containerized instructions.

  ARGUMENTS:
  -h      Flag: Show this help message
  -c      NGSeasy project and run configureation file
  EXAMPLE USAGE:
  aln_novo -c config.file.tsv
EOF 
}

#get options for command line args
  while  getopts "h:c:" opt
  do

      case ${opt} in
	  h)
	  usage #print help
	  exit 0
	  ;;
	  
	  c)
	  config_tsv=${OPTARG}
	  ;;

      esac
  done

#check exists.
  if [[ ! -e ${config_tsv} ]] 
  then
	  echo " ${config_tsv} does not exist "
	  usage;
	  exit 1;
  fi

#run compbio/ngseasy-fastq
sudo docker run \
-d \
-P \
--name aln_novo \
--volumes-from volumes_container \
-t compbio/ngseasy-novoalign:v1.0 /sbin/my_init -- /bin/bash /home/pipeman/ngseasy_scripts/run_ngseasy_novoalign.sh /home/pipeman/ngs_projects/${config_tsv}

}


##--------------------------------------------------##
## NGSeasy: aln_bowtie2
##--------------------------------------------------##

aln_bowtie2 () {

#usage printing func
  usage()
  {
  cat << EOF
  This script sets up the NGSeasy docker Bowtie2 container:
  See NGSEasy containerized instructions.

  ARGUMENTS:
  -h      Flag: Show this help message
  -c      NGSeasy project and run configureation file
  EXAMPLE USAGE:
  aln_bowtie2 -c config.file.tsv
EOF 
}

#get options for command line args
  while  getopts "h:c:" opt
  do

      case ${opt} in
	  h)
	  usage #print help
	  exit 0
	  ;;
	  
	  c)
	  config_tsv=${OPTARG}
	  ;;

      esac
  done

#check exists.
  if [[ ! -e ${config_tsv} ]] 
  then
	  echo " ${config_tsv} does not exist "
	  usage;
	  exit 1;
  fi

#run compbio/ngseasy-fastq
sudo docker run \
-d \
-P \
--name aln_bowtie2 \
--volumes-from volumes_container \
-t compbio/ngseasy-bowtie2:v1.0 /sbin/my_init -- /bin/bash /home/pipeman/ngseasy_scripts/run_ngseasy_bowtie2.sh /home/pipeman/ngs_projects/${config_tsv}

}


##--------------------------------------------------##
## NGSeasy: aln_bwa
##--------------------------------------------------##

aln_bwa () {

#usage printing func
  usage()
  {
  cat << EOF
  This script sets up the NGSeasy docker BWA container:
  See NGSEasy containerized instructions.

  ARGUMENTS:
  -h      Flag: Show this help message
  -c      NGSeasy project and run configureation file
  EXAMPLE USAGE:
  aln_bowtie2 -c config.file.tsv
EOF 
}

#get options for command line args
  while  getopts "h:c:" opt
  do

      case ${opt} in
	  h)
	  usage #print help
	  exit 0
	  ;;
	  
	  c)
	  config_tsv=${OPTARG}
	  ;;

      esac
  done

#check exists.
  if [[ ! -e ${config_tsv} ]] 
  then
	  echo " ${config_tsv} does not exist "
	  usage;
	  exit 1;
  fi

#run compbio/ngseasy-fastq
sudo docker run \
-d \
-P \
--name aln_bwa \
--volumes-from volumes_container \
-t compbio/ngseasy-bwa:v1.0 /sbin/my_init -- /bin/bash /home/pipeman/ngseasy_scripts/run_ngseasy_bwa.sh /home/pipeman/ngs_projects/${config_tsv}

}

##--------------------------------------------------##
## NGSeasy: aln_stampy
##--------------------------------------------------##

aln_stampy () {

#usage printing func
  usage()
  {
  cat << EOF
  This script sets up the NGSeasy docker BWA container:
  See NGSEasy containerized instructions.

  ARGUMENTS:
  -h      Flag: Show this help message
  -c      NGSeasy project and run configureation file
  EXAMPLE USAGE:
  aln_bowtie2 -c config.file.tsv
EOF 
}

#get options for command line args
  while  getopts "h:c:" opt
  do

      case ${opt} in
	  h)
	  usage #print help
	  exit 0
	  ;;
	  
	  c)
	  config_tsv=${OPTARG}
	  ;;

      esac
  done

#check exists.
  if [[ ! -e ${config_tsv} ]] 
  then
	  echo " ${config_tsv} does not exist "
	  usage;
	  exit 1;
  fi

#run compbio/ngseasy-fastq
sudo docker run \
-d \
-P \
--name aln_stampy \
--volumes-from volumes_container \
-t compbio/ngseasy-stampy:v1.0 /sbin/my_init -- /bin/bash /home/pipeman/ngseasy_scripts/run_ngseasy_stampy.sh /home/pipeman/ngs_projects/${config_tsv}

}

  
#### NOTES ##########################################################################################################
## These are the pipeline steps/functions to call

# need a routines to chekc input and output of each step
# if IN exits then DO FUN else exit
# if OUTPUT from last step exists then do next step else run previous step or exit



# at each stage check for sample and project dir and input required
# check if output alread exists and exit if it does then move on 

fastqc_pre
fastq_trimm
fastq_post
aln_novo
aln_stampy
aln_bwa
aln_bowtie2
#aln_gem
#aln_gsnap
#aln_mrsfast
#aln_mrfast
add_read_group
mark_dupes
gatk_clean
gatk_realn
gatk_recal
alignment_qc
var_call
var_call_freebayes
var_call_platypus
var_call_gath_hc
var_call_ensembl
sv_call
sv_call_delly
sv_call_lumpy
sv_call_ensembl
var_anno
var_anno_annovar
var_anno_snpEff
var_anno_ensembl
ngs_reports

  
## TO DO:-

run_ea-ngs.sh at each step is s bash script to perfom NGSeasy step 
run_ea-ngs.sh needs to be copied to a number of worker ngseasy_scripts
eg
run_ea-ngs.sh > run_ngseay_fastqc_pre_aln.sh
run_ea-ngs.sh > run_ngseay_trimmomatic.sh
run_ea-ngs.sh > run_ngseay_fastqc_post_trimmomatic.sh


