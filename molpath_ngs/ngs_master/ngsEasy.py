#!/usr/bin/env python

__doc__=="""
#############################################################################################
# NGS pipeline for molecular diagnostics (can be dockerized)                                #
# -- Organisation: KCL/SLaM/NHS/Viapath                                                     #
# -- Date: 07/08/2014                                                                       #
#############################################################################################
"""
__name__ = 'ngsEasy'
__author__ = "David Brawand, Stephen Newhouse, Amos Folarin, Aditi Gulati"
__credits__ = ['Stephen Newhouse', 'Amos Folarin', 'Aditi Gulati']
__license__ = "LGPL"
__version__ = "0.9"
__maintainer__ = "David Brawand"
__email__ = "dbrawand@nhs.net, stephen.j.newhouse@gmail.com, amosfolarin@gmail.com, aditi.gulati@nhs.net"
__status__ = "Development"  # ["Prototype", "Development",  "Production"]

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   basdic imports and paths
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

import os, sys, re, time
import subprocess

exe_path = os.path.split(os.path.abspath(sys.argv[0]))[0]
#sys.path.insert(0,os.path.abspath(os.path.join(exe_path, "..", "..")))

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   options and input parsing
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
from optparse import OptionParser

usage = "usage: %prog [options] <sample descriptions>"
parser = OptionParser(version="%prog 0.1", usage=usage)

#   configuration files
parser.add_option("-p", dest="pipeconfig", default='pipeline_config.json',metavar="STRING", \
    help="Pipeline parameters (pipeline_config.json)")
parser.add_option("-n", dest="ngsconfig", default='ngs_config.json',metavar="STRING", \
    help="Workflow configuration (ngs_config.json)")

#   general options: verbosity / logging
parser.add_option("-v", dest="verbose", action="count", default=1, \
    help="verbosity (for ruffus meggages)"
     "(eg. -vvv)")
parser.add_option("-l", dest="logfile", default=None, metavar="STRING", \
    help="logfile (default STDERR only)")
parser.add_option("-d" ,dest="debug", action='store_true', default=False,\
    help="display DEBUG messages")

#   pipeline run options
parser.add_option("--touch", dest="touch", default=False, action="store_true", \
    help="just touch")
parser.add_option("--jobs", dest="jobs", default=1, metavar="INT", type="int", \
    help="Specifies the number of jobs (operations) to run in parallel.")
parser.add_option("--flowchart", dest="flowchart", metavar="FILE", type="string", \
                  help="Print flowchart of the pipeline to FILE. Flowchart format "
                       "depends on extension. Alternatives include ('.dot', '.jpg', "
                       "'*.svg', '*.png' etc). Formats other than '.dot' require "
                       "the dot program to be installed (http://www.graphviz.org/).")
parser.add_option("--just_print", dest="just_print", action="store_true", default=False, \
                    help="Only print a trace (description) of the pipeline. "
                         " The level of detail is set by --verbose.")

(options, args) = parser.parse_args()

if not options.flowchart:
    if len(args) == 0:
        parser.print_help()
        parser.error("\nERROR: no sample file provided")


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Logger
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
import logging

# get logger and formatter
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s %(name)-10s %(levelname)-8s %(message)s')

# HANDLERS
# log all to file (evraything!)
if options.logfile:
    logfh = logging.FileHandler(options.logfile)
    logfh.setLevel(logging.DEBUG)
    logfh.setFormatter(formatter)
    logger.addHandler(logfh)
# console logging
logout = logging.StreamHandler(sys.stderr)
if options.debug:
    logout.setLevel(logging.DEBUG)
else:
    logout.setLevel(logging.INFO)
logout.setFormatter(formatter)
logger.addHandler(logout)

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   read requirements/includes
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
import json
from NGSjob import NGSjobs, NGSjob

# read pipeline parameters
logger.info("reading configuration from %s", options.pipeconfig)
with open(options.pipeconfig) as pconf:
    pipeconfig = json.load(pconf)

# read workflows
logger.info("reading workflows from %s", options.ngsconfig)
with open(options.ngsconfig) as nconf:
    ngsconfig = json.load(nconf)

# read sample/ngsjob data
ngsjobs = []
for ngsjobfile in args:
    with open(ngsjobfile) as fh:
        ngsjobs += NGSjobs(fh)

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   setup jobs (collect additional parameters, check files, make directories, write config)
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
inputfiles = []
cfg_file = '.molpath'

# print ngsjob data in readable format
for i, ngsjob in enumerate(ngsjobs):
    targetDir = '/'.join([ pipeconfig["path"]["analysis"], ngsjob.wd() ])
    configpath = targetDir+'/'+cfg_file
    # build configfile and setup directories
    if not os.path.isfile(configpath):
        logger.info('Setting up '+str(ngsjob))
        # get fastq names
        if not ngsjob.FASTQ1:
            try:
                assert ngsjob.projectID and ngsjob.sampleID and ngsjob.runID
            except AssertionError:
                logger.error('Need project, sample and worksheet ID')
                sys.exit(1)
            else:
                ngsjob.FASTQ1 = ngsjob.sampleID + '_' + ngsjob.runID +'_1.fastq'
                ngsjob.FASTQ2 = ngsjob.sampleID + '_' + ngsjob.runID +'_2.fastq'

        # test if FASTQ files exist
        for fastq in ngsjob.fastq(pipeconfig['path']['fastq']):
            try:
                assert os.path.isfile(fastq)
            except:
                logger.error("FASTQ file %s does not exist" % fastq)

        # get platform unit (machine identifier)
        if not ngsjob.RGPU:
            # get ReadGroupPlatformUnit from file and check if reads are paired
            firstlines = []
            for fastq in ngsjob.fastq(pipeconfig['path']['fastq']):
                with open(fastq, 'r') as f:
                    firstlines.append(f.readline().split(' ')[0])
            try:
                assert len(set(firstlines))==1
            except:
                logger.error('Different headers in FastQ files')
                sys.exit(1)

            else:
                #@M00675:4:000000000-A544D:1:1101:19085:2412 1:N:0:3
                ngsjob.RGPU = firstlines[0][1:firstlines[0].find(':')]

        # get runDate
        if not ngsjob.RGDT and ngsjob.sampleSheet:
            try:
                fh = open(ngsjob.sampleSheet)
            except:
                logger.error("cannot find samplesheet file %s" % ngsjob.sampleSheet)
                sys.exit(1)
            else:
                for line in fh:
                    mm = line.match('Date\s+(\d+)\D+(\d+)\D+(\d+)')
                    if mm:
                        if len(mm.group(1))==4 and len(mm.group(3))<=2:
                            ngsjob.RGDT = "-".join(mm.groups())
                        elif len(mm.group(1))<=2 and len(mm.group(3))==4:
                            ngsjob.RGDT = "-".join(mm.groups()[::-1])
                if not ngsjob.RGDT:
                    logger.warn("cannot parse date, use YYYY-MM-DD in samplesheet/RGDT")
            finally:
                fh.close()

        # fallback for RGDT
        if not ngsjob.RGDT:
            logger.warn("no RunDate (RGDT) available, using FastQ file timetamp")
            filetime = lambda x: time.strftime('%Y-%m-%d', time.localtime(min(os.path.getmtime(x),os.path.getctime(x))))
            ngsjob.RGDT = filetime(ngsjobs.fastq(pipeconfig['path']['fastq'])[0])
            print >> sys.stderr, "FASTQ1 last modified: %s" % ngsjob.RGDT  # get earlier date as moving between filesystems could make ctime>mtime

        # analysis check (ngstype is checked after as it can be manually configured)
        try:
            assert ngsjob.ngsAnalysis
            ngsconfig['ngsanalysis'][ngsjob.ngsAnalysis]
        except IndexError:
            logger.error('analysis type unkown')
            sys.exit(1)
        except AssertionError:
            logger.warn('falling back to default analysis parameters')
            ngsjob.ngsAnalysis = 'default'
        except:
            raise

        # ngs data check and set (does not override manually set RG)
        if ngsjob.ngsType:
            try:
                ngsconfig['ngstype'][ngsjob.ngsType]
            except IndexError:
                logger.error("cannot find ngsType %s" % ngsjob.ngsType)
                sys.exit(1)
            except:
                raise
            else:
                for k,v in ngsconfig['ngstype'][ngsjob.ngsType].items():
                    if not getattr(ngsjob, k):
                        setattr(ngsjob, k, v)

        # check if requirements met for pipeline run?
        try:
            missing = ngsjob.unmetRequirements()
            assert len(missing)==0
        except AssertionError:
            logger.error("requirements not met %s" % ','.join(missing))
            sys.exit(1)
        except:
            raise
        else:
            # set RGID if unset
            ngsjob.setRGID()

        # SUCCESS
        #print >> sys.stderr, repr(ngsjob)
        logger.info("Data check successful. Setting up analysis")

        # make target directories
        tmpDir = targetDir + '/tmp'
        try:
            os.makedirs(tmpDir, 0777)
        except OSError as err:
            if err.errno == 17:  # exists
                logger.warn("directory %s exists" % tmpDir)
            else:
                logger.error("Can't create work directory: %s" % err)
                sys.exit(1)

        # symlink with x
        try:
            # first call of inputfiles stores starting files
            symlink = zip(ngsjob.fastq(pipeconfig['path']['fastq']),ngsjob.inputfiles(targetDir))
            for sl in symlink:
                os.symlink(*sl)
        except OSError:
            pass
        except:
            raise

        # write configuration (concludes successful parsing)
        ngsjob.save(configpath)
        logger.info("Configuration complete for %s" % (ngsjob.wd()))

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Start DRMAA session if available
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
try:
    from ruffus.drmaa_wrapper import run_job, error_drmaa_job
    import drmaa
except RuntimeError:  # could not load drmaa
    sge=None
    logger.warn("no DRMAA found -> disabled multiprocessing")
    options.multiprocess = 1

else:
    sge = drmaa.Session()
    sge.initialize()
    logger.info("using DRMAA for job dispatch")


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Functions
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

def loadConfiguration(configpath, configfile=cfg_file):
    # load configuration from file
    configlocation = '/'.join([configpath,configfile])
    try:
        fh = open(configlocation)
    except:
        raise Exception("Could not load configuration (%s)" % (configlocation))
    else:
        ngsjob_config = NGSjob(json.load(fh))
    return ngsjob_config

def run_cmd(cmd_str):
    """
    Throw exception if run command fails
    """
    process = subprocess.Popen(cmd_str, stdout = subprocess.PIPE,
                                stderr = subprocess.PIPE, shell = True)
    stdout_str, stderr_str = process.communicate()
    if process.returncode != 0:
        raise Exception("Failed to run '%s'\n%s%sNon-zero exit status %s" %
                            (cmd_str, stdout_str, stderr_str, process.returncode))

# uses global variables: sge, logger, pipeconfig
def run_sge(cmd_str, jobname, fullwd, cpu=1, mem=2, runlocally=False, scriptdir="tmp" , sge_env={ 'BASH_ENV' : '~/.bashrc' }):
    # get options
    opt = []
    opt.append('-o '+fullwd+'/sge.out')
    opt.append('-e '+fullwd+'/sge.out')
    for k,v in pipeconfig['resources']['gridengine']:
        opt.append(' '.join(k,v))
    if cpu:
        opt.append(' '.join(["-pe multi_thread", str(cpu)]))
    if mem:
        opt.append('-l h_vmem='+str(mem)+'G')
    # dispatch via DRMAA
    try:
        stdout_res, stderr_res  = run_job(
            cmd_str,
            job_name = jobname,
            working_directory = fullwd,
            logger = logger,
            drmaa_session = sge,
            run_locally = runlocally,
            job_other_options = ' '.join(opt),
            job_script_directory = scriptdir,
            job_environment=sge_env, retain_job_scripts = True if scriptdir else False)
    except error_drmaa_job as err:
        raise Exception("\n".join(map(str,["Failed to run:",cmd_str,err,stdout_res,stderr_res])))


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   get input files
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
for i, ngsjob in enumerate(ngsjobs):
    logger.info("Reading configuration from %s" % targetDir)
    targetDir = '/'.join([ pipeconfig["path"]["analysis"], ngsjob.wd() ])
    inputfiles.append(loadConfiguration(targetDir).inputfiles())

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Pipeline tasks
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
from ruffus import *
import inspect  # to get the functions name inside

#--------------------
# FASTQC prefiltering
#--------------------
@transform(inputfiles, formatter("(?P<PAIR>[^\.]+)\.fastq", "(?P<PAIR>[^\.]+)\.fastq"),
            ["{PAIR[0]}_fastqc.html", "{PAIR[1]}_fastqc.html"])
def prefilterFastQC(input_files,output_files):
    # load configurations and name job
    p = loadConfiguration(input_files[0][:input_files[0].rfind('/')])
    cmd = " ".join([ "fastqc", "--noextract", ' '.join(input_files) ])
    # # run on cluster
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=1, mem=2)
    else:
        run_cmd(cmd)


#--------------------
# Adapter trimming
#--------------------
@jobs_limit(4,'iolimit')  # limits jobs with high I/O
@posttask (touch_file( 'trimming.completed' ))
@transform(inputfiles,
            formatter("(?P<PAIR>[^\.]+)\.fastq", "(?P<PAIR>[^\.]+)\.fastq"),
            ["{PAIR[0]}.filtered.fastq", "{PAIR[1]}.filtered.fastq"])
def trimmomatic(input_files, output_files):
    # load configurations
    p = loadConfiguration(input_files[0][:input_files[0].rfind('/')])
    # configure job
    filtered = [ output_files[0], output_files[0].replace('filtered','filtered.unpaired'), \
                 output_files[1], output_files[1].replace('filtered','filtered.unpaired') ]
    cmd = ' '.join([
        pipeconfig['software']['java'],
        '-XX:ParallelGCThreads='+str(pipeconfig['resources']['trimmomatic']['cpu']),
        '-Xmx'+str(pipeconfig['resources']['trimmomatic']['java_mem'])+'g',
        '-jar', pipeconfig['software']['trimmomatic'],
        p.librarytype(),
        '-phred64',
        #'-trimlog', p.RG('SM')+'.trimming.log',
        '-threads', str(pipeconfig['resources']['trimmomatic']['cpu']),
        ' '.join(input_files),
        ' '.join(filtered),
        ' '.join(ngsconfig['ngsanalysis'][p.ngsAnalysis]['trimmomatic']) ])
    # run on cluster
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=1, mem=2)
    # fallback
    else:
        run_cmd(cmd)

#--------------------
# FASTQC postfilter
#--------------------
@transform(trimmomatic, formatter("(?P<PAIR>[^\.]+)\.filtered\.fastq", "(?P<PAIR>[^\.]+)\.filtered\.fastq"),
            ["{PAIR[0]}_fastqc.html", "{PAIR[1]}_fastqc.html"])
def postfilterFastQC(input_files,output_files):
    # load configurations
    p = loadConfiguration(input_files[0][:input_files[0].rfind('/')])
    # create directories
    # qcdir = pipeconfig['path']['analysis']+'/'+p.wd()+"/QC"
    # if not os.path.exists(qcdir):
    #     try: os.mkdir(qcdir)
    #     except: raise
    # cmd = " ".join([ "fastqc", "--noextract", "--outdir="+qcdir, ' '.join(input_files) ])
    cmd = " ".join([ "fastqc", "--noextract", ' '.join(input_files) ])
    # run on cluster
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=1, mem=2)
    # fallback
    else:
        run_cmd(cmd)

#--------------------
# QC report
#--------------------
###### ADD MERGE/COLLATE OF FASTQC and gerneate QC report


#--------------------
# ALIGNMENT
#--------------------
### @collate(animals, regex(r"(.+)\.(.+)\.animal"),  r"\2.results")

@posttask (touch_file( 'alignment.completed' ))
@transform(trimmomatic, formatter("(?P<PAIR>[^\.]+)_1\.filtered\.fastq"), "{PAIR[0]}.sam")
def alignment(input_files, output_file):
    p = loadConfiguration(input_files[0][:input_files[0].rfind('/')])
    # check aligner
    aligner = ngsconfig['ngsanalysis'][p.ngsAnalysis]['aligner']
    try:
        assert aligner in ['bwa','bowtie2','stampy','bowtie','bowtie2','novoalign']
    except:
        logger.warn("aligner %s not available. Falling back to stampy." % aligner)
        aligner = 'stampy'
    # configure job
    if aligner.startswith('bwa'):
        cmd = " ".join([
            pipeconfig['software']['bwa'],
            'mem', '-M',
            '-t '+pipeconfig['resources']['bwa']['cpu'],
            pipeconfig['reference']['genome']['bwaindex'],
            ' '.join(input_files),
            '>', output_file ])
    elif aligner.startswith('stampy'):
        bwa = output_file.replace('.sam','.bwa.sam')
        cmd = " ".join([
            pipeconfig['software']['bwa'],
            'mem', '-M',
            '-t '+pipeconfig['resources']['bwa'],
            pipeconfig['reference']['genome']['sequence'],
            ' '.join(input_files),
            '|', pipeconfig['software']['bwa'], 'view -Sb -', '>', bwa,
            '&&',
            pipeconfig['software']['stampy'],
            '-g', pipeconfig['reference']['genome']['stampyindex'],
            '-h', pipeconfig['reference']['genome']['stampyindex'],
            #'-t', pipeconfig['resources']['stampy']['cpu'],
            #'--bamsortprefix', input_files[0][:input_files[0].rfind('/')]+'/tmp',
            '--bamkeepgoodreads',
            '-M', bwa,
            '-f sam',
            '-o', output_file
            ])
    elif aligner.startswith('novoalign'):
        stat_file = output_file.replace('.sam','.stat')
        cmd = " ".join([
            pipeconfig['software']['novoalign'],
            '-d', pipeconfig['reference']['genome']['novoindex'],
            '-f', ' '.join(input_files),
            '-F STDFQ --Q2Off --3Prime -g 40 -x 6 -r All -i PE 300,150',
            '-c', pipeconfig['resources']['novo']['cpu'],
            '-k -K', stat_file,
            '-o SAM' > output_file
            ])
    elif aligner.startswith('bowtie'):
        cmd = " ".join([
            pipeconfig['software']['bowtie2'],
            '-D 15',
            '-R 2',
            '-N 0',
            '-L 22',
            '-i S,1,1.15',
            '-x', pipeconfig['reference']['genome']['bowtie2index'],
            '-1', input_files[0],
            '-2', input_files[1],
            '-S', output_file
            ])
    # run on cluster or locally
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=1, mem=2)
    else:
        run_cmd(cmd)


#--------------------
# process SAM -> BAM
#--------------------
@transform(alignment, suffix('.sam'), '.bam')
def sam2bam(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    cmd = " ".join([ pipeconfig['software']['samtools'], "view", "-Shb", '-o', output_file, input_file ])
    # run on cluster
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=1, mem=2)
    else:
        run_cmd(cmd)

@transform(sam2bam, suffix('.bam'), '.sorted.bam')
def sortSam(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    cmd = " ".join([ pipeconfig['software']['samtools'],
        "sort", "-f",
        '-@ 2',
        '-m 768M',
        input_file, output_file ])
    # run on cluster
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=1, mem=2)
    else:
        run_cmd(cmd)

@transform(sortSam, suffix('.bam'), '.bam')
def indexBam(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    cmd = " ".join([ pipeconfig['software']['samtools'], "index", input_file ])
    # run on cluster
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=1, mem=2)
    else:
        run_cmd(cmd)

#--------------------
# add/replace ReadGroup
#--------------------
@transform(indexBam, suffix('.bam'), '.addrg.bam')
def addReplaceReadGroups(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    cmd = " ".join([
        pipeconfig['software']['java'],
        '-XX:ParallelGCThreads='+str(pipeconfig['resources']['picard']['cpu']),
        '-Xmx'+str(pipeconfig['resources']['picard']['java_mem'])+'g',
        '-jar', pipeconfig['software']['picard'].rstrip('/')+'/AddOrReplaceReadGroups.jar',
        'TMP_DIR='+'/'.join([ pipeconfig["path"]["analysis"], p.wd(), 'tmp' ]),
        'VALIDATION_STRINGENCY=SILENT',
        'MAX_RECORDS_IN_RAM=100000',
        'CREATE_INDEX=true',
        'SORT_ORDER=coordinate',
        'RGID='+p.RG('ID'),
        'RGLB='+p.RG('LB'),
        'RGPL='+p.RG('PL'),
        'RGPU='+p.RG('PU'),
        'RGSM='+p.RG('SM'),
        'RGDT='+p.RG('DT'),
        'INPUT='+input_file,
        'OUTPUT='+output_file
         ])
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=1, mem=2)
    else:
        run_cmd(cmd)

#--------------------
# mark/remove Duplicates
#--------------------
@transform(addReplaceReadGroups, suffix('.bam'), '.dupemk.bam')
def markDuplicates(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    metrics = output_file.replace('.bam','.stat')
    cmd = " ".join([
        pipeconfig['software']['java'],
        '-XX:ParallelGCThreads='+str(pipeconfig['resources']['picard']['cpu']),
        '-Xmx'+str(pipeconfig['resources']['picard']['java_mem'])+'g',
        '-jar', pipeconfig['software']['picard'].rstrip('/')+'/MarkDuplicates.jar',
        'TMP_DIR='+'/'.join([pipeconfig['path']['analysis'], p.wd(), 'tmp']),
        'VALIDATION_STRINGENCY=SILENT',
        'MAX_RECORDS_IN_RAM=100000',
        'CREATE_INDEX=true',
        "REMOVE_DUPLICATES={rmdup}".format(rmdup="true" if ngsconfig['ngsanalysis'][p.ngsAnalysis]['rmdup'] else "false"),
        'ASSUME_SORTED=true',
        'INPUT='+input_file,
        'OUTPUT='+output_file,
        'METRICS_FILE='+metrics
        ])
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=1, mem=2)
    else:
        run_cmd(cmd)

#--------------------
# GATK indelRealign
#--------------------

@transform(markDuplicates, suffix('.bam'), '.IndelRealigner.intervals')
def realignerTargetCreator(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    cmd = " ".join([
        pipeconfig['software']['java'],
        '-Xmx'+str(pipeconfig['resources']['gatk']['RTC']['java_mem'])+'g',
        '-Djava.io.tmpdir='+'/'.join([pipeconfig['path']['analysis'], p.wd(), 'tmp']),
        '-jar', pipeconfig['software']['gatk'],
        '-T', 'RealignerTargetCreator',
        '-R', pipeconfig['reference']['genome']['sequence'],
        '-nt', str(pipeconfig['resources']['gatk']['RTC']['nt']),
        ' '.join([ '-known '+indl for indl in pipeconfig['reference']['indels'] ]),
        '-I', input_file,
        '-o', output_file
        ])
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=pipeconfig['resources']['gatk']['RTC']['cpu'], mem=pipeconfig['resources']['gatk']['RTC']['mem'])
    else:
        run_cmd(cmd)

@transform([markDuplicates,realignerTargetCreator],
    formatter("(?P<PAIR>[^\.]+)\.bam"),
    "{PAIR[0]}.realn.bam")
def indelRealigner(input_files,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    cmd = " ".join([
        pipeconfig['software']['java'],
        '-Xmx'+str(pipeconfig['resources']['gatk']['IR']['java_mem'])+'g',
        '-Djava.io.tmpdir='+'/'.join([pipeconfig['path']['analysis'], p.wd(), 'tmp']),
        '-jar', pipeconfig['software']['gatk'],
        '-T', 'IndelRealigner',
        '-R', pipeconfig['reference']['genome']['sequence'],
        ' '.join([ '-known '+indl for indl in pipeconfig['reference']['indels'] ]),
        '-I', input_files[0],
        '-targetIntervals', input_files[1],
        '-o', output_file,
        '--consensusDeterminationModel', 'USE_READS',
        '-LOD', '0.4'
        ])
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=pipeconfig['resources']['gatk']['IR']['cpu'], mem=ipeconfig['resources']['gatk']['IR']['mem'])
    else:
        run_cmd(cmd)

@transform(indelRealigner, suffix('.bam'), '.bam')
def indexRealignedBam(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    cmd = " ".join([ pipeconfig['software']['samtools'], "index", input_file ])
    # run on cluster
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=1, mem=2)
    else:
        run_cmd(cmd)

#--------------------
# GATK BaseRecalib
#--------------------

@transform(indexRealignedBam, suffix('.bam'), '.firstBaseRecalibrator.table')
def firstBaseRecalibrator(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    cmd = " ".join([
        pipeconfig['software']['java'],
        '-Xmx'+str(pipeconfig['resources']['gatk']['BR']['java_mem'])+'g',
        '-Djava.io.tmpdir='+'/'.join([pipeconfig['path']['analysis'], p.wd(), 'tmp']),
        '-jar', pipeconfig['software']['gatk'],
        '-T', 'BaseRecalibrator',
        '-R', pipeconfig['reference']['genome']['sequence'],
        '-nct', str(pipeconfig['resources']['gatk']['BR']['nct']),
        ' '.join([ '-knownSites '+indl for indl in pipeconfig['reference']['snps'] ]),
        '-I', input_file,
        '-o', output_file
        ])
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=pipeconfig['resources']['gatk']['BR']['cpu'], mem=pipeconfig['resources']['gatk']['BR']['mem'])
    else:
        run_cmd(cmd)

@transform([indexRealignedBam,firstBaseRecalibrator],
    formatter("(?P<PAIR>[^\.]+)\.bam"),
    "{PAIR[0]}.recal.bam")
def printReads(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    cmd = " ".join([
        pipeconfig['software']['java'],
        '-Xmx'+str(pipeconfig['resources']['gatk']['PR']['java_mem'])+'g',
        '-Djava.io.tmpdir='+'/'.join([pipeconfig['path']['analysis'], p.wd(), 'tmp']),
        '-jar', pipeconfig['software']['gatk'],
        '-T', 'PrintReads',
        '-R', pipeconfig['reference']['genome']['sequence'],
        '-nct', str(pipeconfig['resources']['gatk']['PR']['nct']),
        '--baq', 'RECALCULATE',
        '--baqGapOpenPenalty','40',
        '--BQSR', input_files[1],
        '-I', input_file[0],
        '-o', output_file
        ])
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=pipeconfig['resources']['gatk']['PR']['cpu'], mem=pipeconfig['resources']['gatk']['PR']['mem'])
    else:
        run_cmd(cmd)

# index bam file
@transform(printReads, suffix('.bam'), '.bam')
def indexRecalibratedBam(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    cmd = " ".join([ pipeconfig['software']['samtools'], "index", input_file ])
    # run on cluster
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=1, mem=2)
    else:
        run_cmd(cmd)

@transform(indexRecalibratedBam, suffix('.bam'), '.finalBaseRecalibrator.table')
def finalBaseRecalibrator(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    cmd = " ".join([
        pipeconfig['software']['java'],
        '-Xmx'+str(pipeconfig['resources']['gatk']['BR']['java_mem'])+'g',
        '-Djava.io.tmpdir='+'/'.join([pipeconfig['path']['analysis'], p.wd(), 'tmp']),
        '-jar', pipeconfig['software']['gatk'],
        '-T', 'BaseRecalibrator',
        '-R', pipeconfig['reference']['genome']['sequence'],
        '-nct', str(pipeconfig['resources']['gatk']['BR']['nct']),
        ' '.join([ '-knownSites '+indl for indl in pipeconfig['reference']['snps'] ]),
        '-I', input_file,
        '-o', output_file
        ])
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=pipeconfig['resources']['gatk']['BR']['cpu'], mem=pipeconfig['resources']['gatk']['BR']['mem'])
    else:
        run_cmd(cmd)

@transform([firstBaseRecalibrator,finalBaseRecalibrator],
    formatter("(?P<PAIR>[^\.]+)\.firstBaseRecalibrator\.table"),
    ["{PAIR[0]}.AnalyzeCovariates.csv","{PAIR[0]}.analyzeCovariates.pdf"])
def analyzeCovariates(input_files,output_files):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    cmd = " ".join([
        pipeconfig['software']['java'],
        '-Xmx'+str(pipeconfig['resources']['gatk']['AC']['java_mem'])+'g',
        '-Djava.io.tmpdir='+'/'.join([pipeconfig['path']['analysis'], p.wd(), 'tmp']),
        '-jar', pipeconfig['software']['gatk'],
        '-T', 'AnalyzeCovariates',
        '-R', pipeconfig['reference']['genome']['sequence'],
        '-before', input_files[0],
        '-after', input_files[1],
        '-csv', output_files[0],
        '-plots', output_files[1]
        ])
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=pipeconfig['resources']['gatk']['AC']['cpu'], mem=pipeconfig['resources']['gatk']['AC']['mem'])
    else:
        run_cmd(cmd)

#--------------------
# Haplotype Calling UG
#--------------------

@transform(indexRecalibratedBam, suffix('.bam'), '.UG.vcf')
def unifiedGenotyper(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    cmd = " ".join([
        pipeconfig['software']['java'],
        '-Xmx'+str(pipeconfig['resources']['gatk']['UG']['java_mem'])+'g',
        '-Djava.io.tmpdir='+'/'.join([pipeconfig['path']['analysis'], p.wd(), 'tmp']),
        '-jar', pipeconfig['software']['gatk'],
        '-T', 'UnifiedGenotyper',
        '-R', pipeconfig['reference']['genome']['sequence'],
        '-nct', str(pipeconfig['resources']['gatk']['UG']['nct']),
        '-nt', str(pipeconfig['resources']['gatk']['UG']['nt']),
        '-stand_call_conf', '30',
        '-stand_emit_conf', '10',
        '--dbsnp', pipeconfig['reference']['snps'][0],
        '-dcov', '250',
        '--genotype_likelihoods_model', 'BOTH',
        '--genotyping_mode', 'DISCOVERY',
        '--output_mode', 'EMIT_ALL_CONFIDENT_SITES',
        '--annotation', 'AlleleBalance',
        '--annotation', 'BaseCounts',
        '--annotation', 'BaseQualityRankSumTest',
        '--annotation', 'ChromosomeCounts',
        '--annotation', 'ClippingRankSumTest',
        '--annotation', 'Coverage',
        '--annotation', 'FisherStrand',
        '--annotation', 'GCContent',
        '--annotation', 'HaplotypeScore',
        '--annotation', 'HomopolymerRun',
        '--annotation', 'InbreedingCoeff',
        '--annotation', 'LikelihoodRankSumTest',
        '--annotation', 'LowMQ',
        '--annotation', 'MappingQualityRankSumTest',
        '--annotation', 'MappingQualityZero',
        '--annotation', 'QualByDepth',
        '--annotation', 'RMSMappingQuality',
        '--annotation', 'ReadPosRankSumTest',
        '--annotation', 'SpanningDeletions',
        '--annotation', 'TandemRepeatAnnotator',
        '--annotation', 'VariantType',
        '-I', input_file,
        '-o', output_file
        ])
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=pipeconfig['resources']['gatk']['UG']['cpu'], mem=pipeconfig['resources']['gatk']['UG']['mem'])
    else:
        run_cmd(cmd)

@transform(unifiedGenotyper, suffix('.vcf'), '.filtered.vcf')
def siteFilterUG(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    cmd = " ".join([
        pipeconfig['software']['vcftools'],
        '--vcf',input_file,
        '--recode',
        '--recode-INFO-all',
        '--min-meanDP', ngsconfig['ngsanalysis'][p.ngsAnalysis]['vcf']['minDP'],
        '--minQ', ngsconfig['ngsanalysis'][p.ngsAnalysis]['vcf']['minQ'],
        '--stdout', '>', output_file
        ])
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=1, mem=2)
    else:
        run_cmd(cmd)

#--------------------
# Haplotype Calling HC
#--------------------

@transform(indexRecalibratedBam, suffix('.bam'), '.HC.vcf')
def haplotypeCaller(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    cmd = " ".join([
        pipeconfig['software']['java'],
        '-Xmx'+str(pipeconfig['resources']['gatk']['HC']['java_mem'])+'g',
        '-Djava.io.tmpdir='+'/'.join([pipeconfig['path']['analysis'], p.wd(), 'tmp']),
        '-jar', pipeconfig['software']['gatk'],
        '-T', 'UnifiedGenotyper',
        '-R', pipeconfig['reference']['genome']['sequence'],
        '-nct', str(pipeconfig['resources']['gatk']['HC']['nct']),
        '-nt', str(pipeconfig['resources']['gatk']['HC']['nt']),
        '-stand_call_conf', '30',
        '-stand_emit_conf', '10',
        '--dbsnp', pipeconfig['reference']['snps'][0],
        '-dcov', '250',
        #'--emitRefConfidence', 'GVCF',
        '--genotype_likelihoods_model', 'BOTH',
        '--genotyping_mode', 'DISCOVERY',
        '--output_mode', 'EMIT_ALL_CONFIDENT_SITES',
        '--annotation', 'AlleleBalance',
        '--annotation', 'BaseCounts',
        '--annotation', 'BaseQualityRankSumTest',
        '--annotation', 'ChromosomeCounts',
        '--annotation', 'ClippingRankSumTest',
        '--annotation', 'Coverage',
        '--annotation', 'FisherStrand',
        '--annotation', 'GCContent',
        '--annotation', 'HaplotypeScore',
        '--annotation', 'HomopolymerRun',
        '--annotation', 'InbreedingCoeff',
        '--annotation', 'LikelihoodRankSumTest',
        '--annotation', 'LowMQ',
        '--annotation', 'MappingQualityRankSumTest',
        '--annotation', 'MappingQualityZero',
        '--annotation', 'QualByDepth',
        '--annotation', 'RMSMappingQuality',
        '--annotation', 'ReadPosRankSumTest',
        '--annotation', 'SpanningDeletions',
        '--annotation', 'TandemRepeatAnnotator',
        '--annotation', 'VariantType',
        '-I', input_file,
        '-o', output_file
        ])
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=pipeconfig['resources']['gatk']['HC']['cpu'], mem=pipeconfig['resources']['gatk']['HC']['mem'])
    else:
        run_cmd(cmd)

@transform(haplotypeCaller, suffix('.vcf'), '.filtered.vcf')
def siteFilterHC(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    cmd = " ".join([
        pipeconfig['software']['vcftools'],
        '--vcf',input_file,
        '--recode',
        '--recode-INFO-all',
        '--min-meanDP', ngsconfig['ngsanalysis'][p.ngsAnalysis]['vcf']['minDP'],
        '--minQ', ngsconfig['ngsanalysis'][p.ngsAnalysis]['vcf']['minQ'],
        '--stdout', '>', output_file
        ])
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=1, mem=2)
    else:
        run_cmd(cmd)

#--------------------
# platypus haplotype calling
#--------------------
# only run when stampy is used?
@transform(markDuplicates, suffix('.bam'), '.PL.vcf')
def platypusCallVariants(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    cmd = " ".join([
        pipeconfig['software']['platypus'],
        'callVariants',
        '--bamFiles='+input_file,
        '--refFile='+pipeconfig['reference']['genome']['sequence'],
        '--output='+output_file,
        #'--nCPU='+str(pipeconfig['resources']['platypus']['cpu']),
        #'--minReads='+str(ngsconfig['ngsanalysis'][p.ngsAnalysis]['platypus']['minReads']),
        #'--filterDuplicates='+str(ngsconfig['ngsanalysis'][p.ngsAnalysis]['platypus']['filterDuplicates']),
        #'--trimOverlapping=1',
        ])
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=pipeconfig['resources']['gatk']['HC']['cpu'], mem=pipeconfig['resources']['gatk']['HC']['mem'])
    else:
        run_cmd(cmd)

@transform(platypusCallVariants, suffix('.vcf'), '.filtered.vcf')
def siteFilterPL(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    cmd = " ".join([
        pipeconfig['software']['vcftools'],
        '--vcf',input_file,
        '--recode',
        '--recode-INFO-all',
        '--min-meanDP', ngsconfig['ngsanalysis'][p.ngsAnalysis]['vcf']['minDP'],
        '--minQ', ngsconfig['ngsanalysis'][p.ngsAnalysis]['vcf']['minQ'],
        '--stdout', '>', output_file
        ])
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=1, mem=2)
    else:
        run_cmd(cmd)


'''
# vcfmerge
# target metrics (picard)
# 17. BedTools_DepthOfCoverage
# 18. CollectMultipleMetrics
# 19. Table_Annovar
# 20. Cleanup (structures data)

##---------------------- ALIGNMENT QC ----------------------##

# CollectMultipleMetrics
java -XX:ParallelGCThreads=${NCPU} -Xmx6g -jar /usr/local/pipeline/picardtools/picard-tools-1.115/CollectMultipleMetrics.jar \
TMP_DIR=${SOUT}/tmp \
VALIDATION_STRINGENCY=SILENT \
MAX_RECORDS_IN_RAM=100000 \
INPUT=${SOUT}/alignments/${BAM_PREFIX}.bam \
OUTPUT=${SOUT}/reports/${BAM_PREFIX} \
REFERENCE_SEQUENCE=${REFGenomes}/human_g1k_v37.fasta \
PROGRAM=CollectAlignmentSummaryMetrics \
PROGRAM=CollectInsertSizeMetrics \
PROGRAM=QualityScoreDistribution \
PROGRAM=MeanQualityByCycle;

# CollectAlignmentSummaryMetrics
java -XX:ParallelGCThreads=${NCPU} -Xmx6g -jar /usr/local/pipeline/picardtools/picard-tools-1.115/CollectAlignmentSummaryMetrics.jar \
TMP_DIR=${SOUT}/tmp \
VALIDATION_STRINGENCY=SILENT \
MAX_RECORDS_IN_RAM=100000 \
INPUT=${SOUT}/alignments/${BAM_PREFIX}.bam \
OUTPUT=${SOUT}/reports/${BAM_PREFIX}.alignment_summary_metrics_alt \
REFERENCE_SEQUENCE=${REFGenomes}/human_g1k_v37.fasta \
ASSUME_SORTED=true \
METRIC_ACCUMULATION_LEVEL=SAMPLE;

# CollectWgsMetrics
java -XX:ParallelGCThreads=${NCPU} -Xmx6g -jar /usr/local/pipeline/picardtools/picard-tools-1.115/CollectWgsMetrics.jar \
TMP_DIR=${SOUT}/tmp \
VALIDATION_STRINGENCY=SILENT \
MAX_RECORDS_IN_RAM=100000 \
INPUT=${SOUT}/alignments/${BAM_PREFIX}.bam \
OUTPUT=${SOUT}/reports/${BAM_PREFIX}.wgs_coverage \
REFERENCE_SEQUENCE=${REFGenomes}/human_g1k_v37.fasta \
MINIMUM_MAPPING_QUALITY=20 \
MINIMUM_BASE_QUALITY=20 \
COVERAGE_CAP=1000;

awk 'NR>9' ${SOUT}/reports/${BAM_PREFIX}.wgs_coverage > ${SOUT}/reports/${BAM_PREFIX}.wgs_coverage_hist
sed '8q' ${SOUT}/reports/${BAM_PREFIX}.wgs_coverage > ${SOUT}/reports/${BAM_PREFIX}.wgs_coverage_stats

# FindCoveredIntervals
java -Xmx6g -Djava.io.tmpdir=${SOUT}/tmp -jar /usr/local/pipeline/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T FindCoveredIntervals -R ${REFGenomes}/human_g1k_v37.fasta \
-I ${SOUT}/alignments/${BAM_PREFIX}.bam  \
-o ${SOUT}/reports/${BAM_PREFIX}.CoveredIntervals.list \
--coverage_threshold 4;

# FlagStat
java -Xmx6g -Djava.io.tmpdir=${SOUT}/tmp -jar /usr/local/pipeline/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T FlagStat -R ${REFGenomes}/human_g1k_v37.fasta \
-I ${SOUT}/alignments/${BAM_PREFIX}.bam  \
-o ${SOUT}/reports/${BAM_PREFIX}.FlagStat;

# BAM to BED
/usr/local/pipeline/samtools-0.1.19/samtools view -b -h -q 20 -F 1796  ${SOUT}/alignments/${BAM_PREFIX}.bam  | /usr/local/pipeline/bedtools2/bin/bedtools bamtobed  -i stdin > ${SOUT}/reports/${BAM_PREFIX}.bed;

'''


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Print list of tasks
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
if options.just_print:
    pipeline_printout(sys.stdout, [], verbose=options.verbose)


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Print flowchart
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
elif options.flowchart:
    # use file extension for output format
    output_format = os.path.splitext(options.flowchart)[1][1:]
    pipeline_printout_graph (open(options.flowchart, "w"), output_format, [], no_key_legend = True)

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Run Pipeline
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
else:
    #pipeline_run([trimmomatic], forcedtorun_tasks = [trimmomatic], multiprocess = options.jobs, logger = logger, verbose=options.verbose)
    pipeline_run([], forcedtorun_tasks = [], multiprocess = options.jobs, logger = logger, verbose=options.verbose, touch_files_only=options.touch)

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Cleanly end drmaa session
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
if sge:
    sge.exit()


