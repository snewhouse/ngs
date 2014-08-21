#!/usr/bin/env python

from ngsEasy.helpers import SGE


__doc__=='''
shell script based pipeline (adapted from R/bash scripts)
Classic jobs submission (all at once)
'''

def trimAdapters(p,sge):
    jobname = "trimming_"+p.uniqueID()
    runscript = "ngs_AdapterTrimming.sh"
    # send job
    sge.qsub(exe=p.pipeline_dir+'/'+runscript+" "+p.uniqueID()+" "+p.wd()+" "+p.fastq_names(), name=jobname)
    # add dependency for next job
    sge.depends.append(jobname)

def mapReadsNovoalign(p,sge):
    # novoalign
    jobname = "novoalign_"+p.uniqueID()
    runscript = "ngs_novoalign.sh"
    sge.qsub(exe=p.pipeline_dir+'/'+runscript+" "+p.uniqueID()+" "+p.wd(), name=jobname)
    # add dependency for next job
    sge.depends.append(jobname)
    return

def mapReadsStampy(p,sge):
    # novoalign
    jobname = "stampy_"+p.uniqueID()
    runscript = "ngs_stampy.sh"
    sge.qsub(exe=p.pipeline_dir+'/'+runscript+" "+p.uniqueID()+" "+p.wd(), name=jobname)
    # add dependency for next job
    sge.depends.append(jobname)
    return

def sam2bam(p,sge):
    jobname = 'sam2bam.'+p.uniqueID()
    runscript = 'ngs_sam2bam.sh'
    sge.qsub(exe=p.pipeline_dir+'/'+runscript+" "+p.uniqueID()+" "+p.wd(), name=jobname)
    sge.depends.append(jobname)
    return

def sortSam(p, sge):
    jobname = 'sortSam.'+p.uniqueID()
    runscript = 'ngs_sam2bam.sh'
    sge.qsub(exe=p.pipeline_dir+'/'+runscript+" "+p.uniqueID()+" "+p.wd()+" "+p.tmp(), name=jobname)
    sge.depends.append(jobname)
    return

def addReplaceReadGroups(p, sge):
    jobname = 'AddOrReplaceReadGroups.'+p.uniqueID()
    runscript = 'ngs_AddOrReplaceReadGroups.sh'
    sge.qsub(exe=p.pipeline_dir+'/'+runscript+" "+\
        " ".join([p.uniqueID(), p.wd(), p.tmp(), p.RGID(), p.RGLB(), p.RGPL(), p.RGPU(), p.RGSM(), p.RGCN(), p.RGDS(), p.RGDT()]), \
        name=jobname)
    sge.depends.append(jobname)
    return

def MarkDuplicates(p, sge):
    jobname = 'AddOrReplaceReadGroups.'+p.uniqueID()
    runscript = 'ngs_AddOrReplaceReadGroups.sh'
    sge.qsub(exe=p.pipeline_dir+'/'+runscript+" "+\
        " ".join([p.uniqueID(), p.wd(), p.tmp(), p.RGID(), p.RGLB(), p.RGPL(), p.RGPU(), p.RGSM(), p.RGCN(), p.RGDS(), p.RGDT()]), \
        name=jobname)
    sge.depends.append(jobname)
    return


#----------------------------------------------------------------------#
# 5. MarkDuplicates
#----------------------------------------------------------------------#

###echo ">>>>>" `date` " :-> " "Running MarkDuplicates"

qsub -o ${SGE_OUT} -e ${SGE_OUT} -q ${queue_name} -N MarkDuplicates.${sample_name} -hold_jid AddOrReplaceReadGroups.${sample_name} -l h_vmem=${sge_h_vmem}G  -pe multi_thread 1 -M ${email_contact} -m beas ${ngs_pipeline}/ngs_MarkDuplicates.sh ${sample_name} ${sample_dir} ${sample_temp};

##############################
## END ALIGNMENT PIPELINE   ##
##############################

#########################
## BEGIN GATK CLEANING ##
#########################

#----------------------------------------------------------------------#
# 7. RealignerTargetCreator
#----------------------------------------------------------------------#

##echo ">>>>>" `date` " :-> " "Running RealignerTargetCreator"

qsub -o ${SGE_OUT} -e ${SGE_OUT} -q ${queue_name} -N RealignerTargetCreator.${sample_name} -hold_jid MarkDuplicates.${sample_name} -l h_vmem=${gatk_h_vmem}G -M ${email_contact} -m beas ${ngs_pipeline}/ngs_RealignerTargetCreator.sh ${sample_name} ${sample_dir} ${sample_temp};


#----------------------------------------------------------------------#
# 8. IndelRealigner
#----------------------------------------------------------------------#

##echo ">>>>>" `date` " :-> " "Running IndelRealigner"

qsub -o ${SGE_OUT} -e ${SGE_OUT} -q ${queue_name} -N IndelRealigner.${sample_name} -hold_jid RealignerTargetCreator.${sample_name} -l h_vmem=${gatk_h_vmem}G -M ${email_contact} -m beas ${ngs_pipeline}/ngs_IndelRealigner.sh ${sample_name} ${sample_dir} ${sample_temp};


#----------------------------------------------------------------------#
# 9. BaseRecalibrator before recal
#----------------------------------------------------------------------#

##echo ">>>>>" `date` " :-> " "Running BaseRecalibrator before QUAL SCORE RECALIBRATION"

qsub -o ${SGE_OUT} -e ${SGE_OUT} -q ${queue_name} -N BaseRecalibrator_before.${sample_name} -hold_jid IndelRealigner.${sample_name} -l h_vmem=${gatk_h_vmem}G -M ${email_contact} -m beas ${ngs_pipeline}/ngs_BaseRecalibrator_before.sh ${sample_name} ${sample_dir} ${sample_temp};


#----------------------------------------------------------------------#
# 10. PrintReads = QUAL SCORE RECALIBRATION
#----------------------------------------------------------------------#

###echo ">>>>>" `date` " :-> " "Running PrintReads > QUAL SCORE RECALIBRATION "

qsub -o ${SGE_OUT} -e ${SGE_OUT} -q ${queue_name} -N PrintReads_BQSR.${sample_name} -hold_jid BaseRecalibrator_before.${sample_name} -l h_vmem=${gatk_h_vmem}G -M ${email_contact} -m beas ${ngs_pipeline}/ngs_PrintReads_BQSR.sh ${sample_name} ${sample_dir} ${sample_temp};


#----------------------------------------------------------------------#
# 11. BaseRecalibrator after recal
#----------------------------------------------------------------------#

##echo ">>>>>" `date` " :-> " "Running BaseRecalibrator after QUAL SCORE RECALIBRATION"

############ qsub -o ${SGE_OUT}  -e ${SGE_OUT}  -q ${queue_name} -N BaseRecalibrator_after.${sample_name} -hold_jid PrintReads_BQSR.${sample_name} -l h_vmem=${gatk_h_vmem}G -M ${email_contact} -m beas ${ngs_pipeline}/ngs_BaseRecalibrator_after.sh ${sample_name} ${sample_dir} ${sample_temp};


#----------------------------------------------------------------------#
# 12. AnalyzeCovariates before & after recal
#----------------------------------------------------------------------#
##echo ">>>>>" `date` " :-> " "Running AnalyzeCovariates"
### qsub -o ${SGE_OUT}  -e ${SGE_OUT}  -q ${queue_name} -N AnalyzeCovariates_before_and_after_BQSR.${sample_name} -hold_jid PrintReads_BQSR.${sample_name} -l h_vmem=${gatk_h_vmem}G -M ${email_contact} -m beas ${ngs_pipeline}/ngs_AnalyzeCovariates_before_and_after_BQSR.sh \
#### ${sample_name} ${sample_dir} ${sample_temp};

##############################################
##  CALL VARIANTS SINGLE SAMPLE ##############
##############################################

#----------------------------------------------------------------------#
# 13. HaplotypeCaller
#----------------------------------------------------------------------#
#<<<<<<< HEAD

#=======
#>>>>>>> fefd264974fa61b46813acfa885d801845d5fffc
qsub -o ${SGE_OUT} -e ${SGE_OUT} -q ${queue_name} -N HaplotypeCaller.${sample_name} -hold_jid PrintReads_BQSR.${sample_name} -l h_vmem=${gatk_h_vmem}G -M ${email_contact} -m beas ${ngs_pipeline}/ngs_HaplotypeCaller.sh \
${sample_name} ${sample_dir} ${sample_temp} 30 10;

#----------------------------------------------------------------------#
# 14. VCFtoolsSiteFilter_on_HaplotypeCaller_Output
#----------------------------------------------------------------------#

qsub -o ${SGE_OUT} -e ${SGE_OUT} -q ${queue_name} -N VCFtoolsSiteFilter_HaplotypeCaller.${sample_name} -hold_jid HaplotypeCaller.${sample_name} -l h_vmem=${gatk_h_vmem}G -M ${email_contact} -m beas ${ngs_pipeline}/ngs_VCFtoolsSiteFilter_HaplotypeCaller.sh \
${sample_name} ${sample_dir} ${sample_temp};

#----------------------------------------------------------------------#
# 15. UnifiedGenotyper
#----------------------------------------------------------------------#

qsub -o ${SGE_OUT} -e ${SGE_OUT} -q ${queue_name} -N UnifiedGenotyper.${sample_name} -hold_jid PrintReads_BQSR.${sample_name} -l h_vmem=${gatk_h_vmem}G -M ${email_contact} -m beas ${ngs_pipeline}/ngs_UnifiedGenotyper.sh \
${sample_name} ${sample_dir} ${sample_temp} 30 10;

## to do : discovery > list > merge novels with knowns and re-genotype. Include a genotype all bases option

#----------------------------------------------------------------------#
# 16.VCFtoolsSiteFilter_on_UnifiedGenotyper_Output
#----------------------------------------------------------------------#

qsub -o ${SGE_OUT} -e ${SGE_OUT} -q ${queue_name} -N VCFtoolsSiteFilter_UnifiedGenotyper.${sample_name} -hold_jid UnifiedGenotyper.${sample_name} \
-l h_vmem=${gatk_h_vmem}G -M ${email_contact} -m beas ${ngs_pipeline}/ngs_VCFtoolsSiteFilter_UnifiedGenotyper.sh ${sample_name} ${sample_dir} ${sample_temp};

#############################################################
## summary metrics ##########################################
#############################################################

#----------------------------------------------------------------------#
# 17. BedTools_DepthOfCoverage
#----------------------------------------------------------------------#

qsub -o ${SGE_OUT} -e ${SGE_OUT} -q ${queue_name} -N DepthOfCoverage.${sample_name} -hold_jid PrintReads_BQSR.${sample_name} -l h_vmem=${gatk_h_vmem}G -M ${email_contact} -m beas ${ngs_pipeline}/ngs_BedTools_DepthOfCoverage.sh \
${sample_name} ${sample_dir} ${sample_temp} ${bed_list} ${bed_type};


#----------------------------------------------------------------------#
# 18. CollectMultipleMetrics
#----------------------------------------------------------------------#

qsub -o ${SGE_OUT} -e ${SGE_OUT} -q ${queue_name} -N CollectMultipleMetrics.${sample_name} -hold_jid PrintReads_BQSR.${sample_name} -l h_vmem=${gatk_h_vmem}G -M ${email_contact} -m beas ${ngs_pipeline}/ngs_CollectMultipleMetrics.sh \
${sample_name} ${sample_dir} ${sample_temp}



#######################################################################
## Variant annotations ################################################
#######################################################################

#----------------------------------------------------------------------#
# 19. Table_Annovar
#----------------------------------------------------------------------#

qsub -o ${SGE_OUT} -e ${SGE_OUT} -q ${queue_name} -N annovar_HaplotypeCaller.${sample_name} -hold_jid VCFtoolsSiteFilter_HaplotypeCaller.${sample_name} \
-l h_vmem=${gatk_h_vmem}G -M ${email_contact} -m beas ${ngs_pipeline}/ngs_table_annovar_HaplotypeCaller_hg19.sh ${sample_name} ${sample_dir} ${sample_temp} "HaplotypeCaller";

qsub -o ${SGE_OUT} -e ${SGE_OUT} -q ${queue_name} -N annovar_UnifiedGenotyper.${sample_name} -hold_jid VCFtoolsSiteFilter_UnifiedGenotyper.${sample_name} \
-l h_vmem=${gatk_h_vmem}G -M ${email_contact} -m beas ${ngs_pipeline}/ngs_table_annovar_UnifiedGenotyper_hg19.sh ${sample_name} ${sample_dir} ${sample_temp} "UnifiedGenotyper";

#########################
## END ##################
#########################




##########################################################################################
##### OLD STUFF ##########################################################################
##########################################################################################
# read mapping (PE/SE, illumina/iontorrent)
def mapReadsLegacy(p,sge):
    jobname = "novoalign_"+p.uniqueID()
    # determine script to run
    if p.RGPL().startswith('illumina'):
        if p.isPE():
            logging("Reads are PE ILLUMINA")
            runscript = "ngs_novoalign.illumina."+p.QUAL+".PE.sh"
        else:
            logging("Reads are SE ILLUMINA")
            runscript = "ngs_novoalign.illumina."+p.QUAL+".SE.sh"

    elif p.RGPL().startswith("iontorrent"):
        try:
            assert p.isPE()
        except AssertionError:
            raise Exception("IONTORRENT cannot be PE")
        else:
            logging("Reads are SE IONTORRENT")
            runscript = "ngs_novoalign.IONTORRENT."+p.QUAL+".SE.sh"
    else:
        raise Exception("mapping configuration????")

    # send job
    sge.qsub(exe=p.pipeline_dir+'/'+runscript+" "+p.fastq_names()+" "+p.uniqueID()+" "+p.wd(), name=jobname)

    # add dependency for next job
    sge.depends.append(jobname)

    return

${ngs_pipeline}/ngs_sam2bam.sh ${sample_name} ${sample_dir};

${ngs_pipeline}/ngs_SortSam.sh ${sample_name} ${sample_dir} ${sample_temp};
