#!/bin/bash
# Collection of NGSeasy Functions
# Stephen Newhouse <stephen.j.newhouse@gmail.com>
# Version 0.9.0

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

# run it
  

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
--name fastqc_and_trim \
--volumes-from volumes_container \
-t compbio/ngseasy-fastqc:v1.0 /sbin/my_init -- /bin/bash /home/pipeman/ngseasy_scripts/run_ea-ngs.sh /home/pipeman/ngs_projects/${config_tsv}

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

sudo docker run \
-d \
-P \
--name fastqc_and_trim \
--volumes-from volumes_container \
-t compbio/ngseasy-fastqc:v1.0 /sbin/my_init -- /bin/bash /home/pipeman/ngseasy_scripts/run_ea-ngs.sh /home/pipeman/ngs_projects/${config_tsv}

}


##--------------------------------------------------##
## NGSeasy: FastQC Post-Trimmomatic
##--------------------------------------------------##




  
#### NOTES ##########################################################################################################
## These are the pipeline steps/functions to call

need a routines to chekc input and output of each step
if IN exits then DO FUN else exit
if OUTPUT from last step exists then do next step else run previous step or exit

fastqc_pre
fastq_trimm
fastq_post
aln_novo
aln_stampy
aln_bwa
aln_bowtie2
aln_gem
aln_gsnap
aln_mrsfast
aln_mrfast
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


