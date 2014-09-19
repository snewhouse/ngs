#!/bin/bash
# NSGeasy Project Initiation
# Makes Project and Sample Directories needed for NGSesy
# Stephen Newhouse <stephen.j.newhouse@gmail.com>
# Amos Folarin <afolarin@gmail.com>


#------------------------------------------------------------#
# Read ngs config file 
#------------------------------------------------------------#
config_file=${1}
host_vol_dir=${2}

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

## Print input 
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

 
#------------------------------------------------------------#
# Checking Project Directory 
#------------------------------------------------------------#
  
echo  "#------------------------------------------------------------#"
echo  " Checking Project/Sample Directory Structure"
echo  "#------------------------------------------------------------#"
echo ""
  
  if [ ! -d $PROJECT_DIR/${POJECT_ID} ]
    then
    echo  " Project Directory Does not exist" 
    echo  " Setting up correct Project Directory"
    echo ""
    echo  " Making Directory : [$PROJECT_DIR/${POJECT_ID}/]"
    echo  " Making Directory : [$PROJECT_DIR/${POJECT_ID}/config_files]"
    echo  " Making Directory : [$PROJECT_DIR/${POJECT_ID}/cohort_vcfs]"
    echo ""
    
        mkdir ${PROJECT_DIR}/${POJECT_ID}/
        mkdir ${PROJECT_DIR}/${POJECT_ID}/config_files
        mkdir ${PROJECT_DIR}/${POJECT_ID}/cohort_vcfs
        
  else
    echo  " Project Directory Exists"
  fi

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
echo ""
DATE=`date +"%d%m%y"`
echo  " Run Date : [$DATE]"
echo  " Project_Id : [$POJECT_ID]" 
echo  " Sample_Id  : [$SAMPLE_ID]"
echo  " Project Directory : [$PROJECT_DIR]"
echo  " Platform : [$NGS_PLATFORM]"
echo  " NGS Type : [$NGS_TYPE]"
echo  " Clean up TRUE/FALE : [$CLEANUP]"
echo  " Number of cpu : [$NCPU]"
echo  " Output Directory : [${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}]"
echo ""

echo  "#------------------------------------------------------------#"
echo  " Call NGS Pipeline "
echo  "#------------------------------------------------------------#"
echo ""
 
echo  " Saving Project configuration settings to : [${PROJECT_DIR}/${POJECT_ID}/config_files/master.config.file]"
echo  "$DATE $f1 $f2 $f3 $f4 $f5 $f6 $f7 $f8 $f9 $f10 $f11 $f12 $f13 $f14 $f15 ${SAMPLE_ID}.${NGS_TYPE}.${NGS_PLATFORM}.${ALIGNER}.${DATE} ${PROJECT_DIR}/${POJECT_ID}" >> ${PROJECT_DIR}/${POJECT_ID}/config_files/project.config.file;
  
done < ${config_file}


