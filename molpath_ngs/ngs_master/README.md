
reference variants have to be in order: DBSNP,HapMap,1000Gomni!!! 




Pipeline Bash Scripts
====================================
- ...coming soon...
- easyngs_full.v1.0.sh
- aln,gatk_clean,var-call,sv-call,var-anno,full,full_nogatk

Descriptions  
config file set up  


molpath.py parameter dependencies

## parameters
### sample file
projectID
sampleID
worksheetID
    RGID
    FASTQ1
    FASTQ2
        RGPU
sampleSheet(FILE)
    RGDT
ngsType
    RGLB
    RGCN
    RGPL
    RGDS
ngsAnalysis
    cleanupIntermediary

# FASTQ naming convention
 ngsjob.FASTQ1 = ngsjob.sampleID + '_' + ngsjob.worksheetID +'_'+'_1.fastq'

# field aliases
runDate: RGDT
runID: worksheetID

BEDtype
WGS TAS EXO

for normal/tumor cancer panels use:
    Platypus.py --assemble=1

Variantfilter:
    bcftools filter


