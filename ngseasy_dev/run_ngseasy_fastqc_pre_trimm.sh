#!/bin/bash

# run_ngseasy_fastqc_pre_trimm.sh 

run_type=${1}
config_tsv=${2}

# usage function 
usage()
  {
  cat << EOF
  This script sets up the NGSeasy docker fastqc container:
  See NGSEasy containerized instructions.

  EXAMPLE USAGE:
  bash run_ngseasy_fastqc_pre_trimm.sh <config.file.tsv>
EOF 
}

if [ "${run_type}" == "config" ]
then

#-------------------------------------------------------------------#

#check exists.
  if [[ ! -e ${config_tsv} ]] 
  then
	  echo " ${config_tsv} does not exist "
	  usage;
	  exit 1;
  fi

#-------------------------------------------------------------------#

# read config file 

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

#-------------------------------------------------------------------#
#BAM PREFIX 
  BAM_PREFIX=${SAMPLE_ID}.${NGS_TYPE}.${NGS_PLATFORM}.${ALIGNER}
  echo " NGSeasy: Setting BAM_PREFIX directory [$BAM_PREFIX]"
 
#OUTPUT SAMPLE DIR
 SOUT=${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}
#-------------------------------------------------------------------#


if [ ! -e ${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID} ]
then
  echo " NGSeasy: Cant Find Project directory. This is then end. Please Stop and check everything is ok " `date`
  exit 1

else 
  echo " NGSeasy: Setting OUTPUT directory [${SOUT}]"
fi

#copy fastq files to sample directory
if [ ! -s ${SOUT}/fastq/${FASTQ1} ] && [ ! -s ${SOUT}/fastq/${FASTQ2} ]
then
  echo " NGSeasy: Copying fastq files from ${FASTQDIR}/ to ${SOUT}/fastq/ " `date`
    cp ${FASTQDIR}/${FASTQ1} ${SOUT}/fastq/${FASTQ1}
    cp ${FASTQDIR}/${FASTQ2} ${SOUT}/fastq/${FASTQ2}

else
  echo " NGSeasy: Fastq Files exist in  ${SOUT}/fastq/ " `date`
  ls ${SOUT}/fastq/
fi

#set new names for copied fastq files
  rawFASTQ1=`basename ${SOUT}/fastq/${FASTQ1} _1.fq.gz`
  rawFASTQ2=`basename ${SOUT}/fastq/${FASTQ2} _2.fq.gz`
    
echo " NGSeasy: Fastq Basename : [$rawFASTQ1] "

#-------------------------------------------------------------------#
#FASTQC on raw files
echo " NGSeasy: START Pre-Alignment QC " `date`

#check if qc'd data alread exists 
if [ ! -s ${SOUT}/fastq/${rawFASTQ1}_1.fq_fastqc.zip ] && [ ! -s ${SOUT}/fastq/${rawFASTQ2}_2.fq_fastqc.zip ]
then
  echo " NGSeasy: Run Pre-Alignment QC on raw Fastq files " `date`

  /usr/local/pipeline/FastQC/fastqc \
    --threads ${NCPU} \
    --extract \
    --quiet \
    --dir ${SOUT}/tmp \
    --outdir ${SOUT}/fastq \
    ${SOUT}/fastq/${rawFASTQ1}_1.fq.gz \
    ${SOUT}/fastq/${rawFASTQ2}_2.fq.gz
    
else
  echo " NGSeasy: Pre-Alignment QC on raw Fastq files Already run " `date`
fi

echo " NGSeasy: END Pre-Alignment QC  " `date`

done < ${config_file}

else
    POJECT_ID=$f1
    SAMPLE_ID=$f2
    FASTQ1=$f3
    FASTQ2=$f4
    PROJECT_DIR=$f5
#OUTPUT SAMPLE DIR
 SOUT=${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}
#-------------------------------------------------------------------#


if [ ! -e ${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID} ]
then
  echo " NGSeasy: Cant Find Project directory. This is then end. Please Stop and check everything is ok " `date`
  exit 1

else 
  echo " NGSeasy: Setting OUTPUT directory [${SOUT}]"
fi

#copy fastq files to sample directory
if [ ! -s ${SOUT}/fastq/${FASTQ1} ] && [ ! -s ${SOUT}/fastq/${FASTQ2} ]
then
  echo " NGSeasy: Copying fastq files from ${FASTQDIR}/ to ${SOUT}/fastq/ " `date`
    cp ${FASTQDIR}/${FASTQ1} ${SOUT}/fastq/${FASTQ1}
    cp ${FASTQDIR}/${FASTQ2} ${SOUT}/fastq/${FASTQ2}

else
  echo " NGSeasy: Fastq Files exist in  ${SOUT}/fastq/ " `date`
  ls ${SOUT}/fastq/
fi

#set new names for copied fastq files
  rawFASTQ1=`basename ${SOUT}/fastq/${FASTQ1} _1.fq.gz`
  rawFASTQ2=`basename ${SOUT}/fastq/${FASTQ2} _2.fq.gz`
    
echo " NGSeasy: Fastq Basename : [$rawFASTQ1] "

#-------------------------------------------------------------------#
#FASTQC on raw files
echo " NGSeasy: START Pre-Alignment QC " `date`

#check if qc'd data alread exists 
if [ ! -s ${SOUT}/fastq/${rawFASTQ1}_1.fq_fastqc.zip ] && [ ! -s ${SOUT}/fastq/${rawFASTQ2}_2.fq_fastqc.zip ]
then
  echo " NGSeasy: Run Pre-Alignment QC on raw Fastq files " `date`

  /usr/local/pipeline/FastQC/fastqc \
    --threads ${NCPU} \
    --extract \
    --quiet \
    --dir ${SOUT}/tmp \
    --outdir ${SOUT}/fastq \
    ${SOUT}/fastq/${rawFASTQ1}_1.fq.gz \
    ${SOUT}/fastq/${rawFASTQ2}_2.fq.gz
    
else
  echo " NGSeasy: Pre-Alignment QC on raw Fastq files Already run " `date`
fi

echo " NGSeasy: END Pre-Alignment QC  " `date`

fi


