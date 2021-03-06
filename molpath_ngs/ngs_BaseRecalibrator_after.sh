#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -p -0.99999999999999999999999999999999999999999999999999
#$ -V

###############################################
## BaseRecalibrator with default covars
################################################

sample_name=${1}
sample_dir=${2}
sample_temp=${3}

cd ${sample_dir}

${java_v1_7}/java  -Xmx${gatk_java_mem}g -Djava.io.tmpdir=${sample_temp} -jar ${ngs_gatk}/GenomeAnalysisTK.jar -T BaseRecalibrator -R ${reference_genome_seq}  \
-I ${sample_dir}/${sample_name}.novorecal.bam \
-o ${sample_dir}/${sample_name}.novorecal.recal_data.table \
-knownSites ${b37_1000G_indels} \
-knownSites ${b37_Mills_Devine_2hit_indels} \
-knownSites ${b37_1000G_omni2_5} \
-knownSites ${b37_1000G_snps} \
-knownSites ${b37_hapmap_3_3} \
-knownSites ${b37_dbsnp};

##########mv -v BaseRecalibrator_after.${sample_name}.* ${sample_dir}/sge_out/


#-nct 1 \
#-log ${sample_dir}/${sample_name}.novorealn.BaseRecalibrator.log;

#--covariate / -cov ( String[] )
#One or more covariates to be used in the recalibration. 
#Can be specified multiple times. 
#Note that the ReadGroup and QualityScore covariates are required and do not need to be specified. 
##Also, unless --no_standard_covs is #specified, the Cycle and Context covariates are standard and are included by default.
#Use the --list argument to see the available covariates.
