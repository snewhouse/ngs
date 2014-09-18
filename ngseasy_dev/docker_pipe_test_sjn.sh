#!/bin/bash

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

### fastqc YES or No

fastqc_pre () {
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






### run it
  ${fastqc_and_trim}    

### trimming YES or No  
  fastq_trimm=$(sudo docker run \
		    -d \
		    -P \
		    --name fastqc_and_trim \
		    --volumes-from volumes_container \
		    -t compbio/ngseasy-fastqc:v1.0 /sbin/my_init -- /bin/bash /home/pipeman/ngseasy_scripts/run_ea-ngs.sh /home/pipeman/ngs_projects/${config_tsv}
		    )
  
### run it
  ${fastq_trimm}  
  
  
#### NOTES
## These are the pipeline steps

fastqc_pre
fastq_trimm
fastq_post
alignment
add_read_group
mark_dupes
gatk_clean
alignment_qc
var_call
sv_call
var_anno
ngs_reports

  
  
### Alignment 
  aln_novo=$(sudo docker  run -i -t compbio/)
  aln_stampy=$(sudo docker  run -i -t compbio/)
  aln_bwa=$(sudo docker  run -i -t compbio/)
  aln_bowtie2=$(sudo docker  run -i -t compbio/)
  aln_gem=$(sudo docker  run -i -t compbio/)
  aln_gsnap=$(sudo docker  run -i -t compbio/)
  aln_mrsfast=$(sudo docker  run -i -t compbio/)
  aln_mrfast=$(sudo docker  run -i -t compbio/)


## Add read groups
  add_read_group=$(sudo docker  run -i -t compbio/)

## MarkDuplicates
  mark_dupes=$(sudo docker  run -i -t compbio/)


### GATK Clean YES or NO
  gatk_clean=$(sudo docker  run -i -t compbio/)

### Alignment QC
  alignment_qc=$(sudo docker  run -i -t compbio/)


### SNV Calling
  var_call=$(sudo docker  run -i -t compbio/) # call all 3 haplotype based callers - GATK, Freebayes and Platypus - consenus as truth for GATK variant recalibration routine

### SV Calling
  sv_call=$(sudo docker  run -i -t compbio/)

### Variant Annotation
  var_anno=$(sudo docker  run -i -t compbio/)

### Report Summaries 


##--------------------------------------------------##
## END ALL CONTAINERS 
##--------------------------------------------------##
  sudo docker kill ${volumes_container}


