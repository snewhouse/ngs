#!/usr/bin/env python

__doc__=="""
#############################################################################################
# NGS pipeline for molecular diagnostics (can be dockerized)                                #
# -- Organisation: KCL/SLaM/NHS/Viapath                                                     #
# -- Date: 07/08/2014                                                                       #
#############################################################################################
"""
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
parser.add_option("-v", dest="verbose", action="count", default=0, \
    help="Print more detailed messages (eg. -vvv)")
parser.add_option("-l", dest="logfile", default=None, metavar="STRING", \
    help="logfile (default STDERR)")

#   pipeline run options
parser.add_option("--jobs", dest="jobs", default=2, metavar="INT", type="int", \
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
logging.basicConfig(
            filename=options.logfile,
            filemode='w',
            format='%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%m/%d/%Y %I:%M:%S %p',
            level=logging.DEBUG)
logger = logging.getLogger("MolPath")

if options.verbose:
    stderrhandler = logging.StreamHandler(sys.stderr)
    stderrhandler.setFormatter(logging.Formatter("    %(message)s"))
    stderrhandler.setLevel(logging.DEBUG)
    logger.addHandler(stderrhandler)


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

# print ngsjob data in readable format
for i, ngsjob in enumerate(ngsjobs):
    print i, ngsjob

    # get fastq names
    try:
        assert self.projectID and sampleID and worksheetID
    except AssertionError:
        logger.error('Need project, sample and worksheet ID')
    else:
        if not ngsjob.FASTQ1:
            ngsjob.FASTQ1 = ngsjob.sampleID + '_' + ngsjob.worksheetID +'_1.fastq'
            ngsjob.FASTQ2 = ngsjob.sampleID + '_' + ngsjob.worksheetID +'_2.fastq'

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
        else:
            #@M00675:4:000000000-A544D:1:1101:19085:2412 1:N:0:3
            ngsjob.RGPU = firstlines[0][1:firstlines[0].find(':')]

    if not ngsjob.RGDT and ngsjob.sampleSheet:
        try:
            fh = open(ngsjob.sampleSheet)
        except:
            logger.error("cannot find samplesheet file %s" % ngsjob.sampleSheet)
        else:
            for line in fh:
                mm = line.match('Date\s+(\d+)\D+(\d+)\D+(\d+)')
                if mm:
                    if len(mm.group(1))==4 and len(mm.group(3))<=2:
                        ngsjob.RGDT = mm.group(1) +
                    elif:
                        ngsjob.RGDT = mm.group(1)
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

    # check if requirements met for pipeline run?
    try:
        missing = ngsjob.unmetRequirements()
        assert len(missing)==0
    except AssertionError:
        logger.error("requirements not met %s" % ','.join(missing))

    logger.info("Data check successful. Setting up analysis")

    # make target directories and write configuration
    targetDir = '/'.join([ pipeconfig["paths"]["analysis"], ngsjob.wd() ])
    try:
        os.makedirs(targetDir, 0777)
    except OSError as err:
        if "[Errno 17]" in err:  # exists
            logger.warn("directory %s exists" % targetDir)
        else:
            logger.error("directory %s creation error: %s" % (targetDir, err))
    else:
        ngsjob.save(targetDir+'/.molpath')

    # symlink with x
    try:
        symlink = zip(ngsjob.fastq(pipeconfig['path']['fastq']),ngsjob.fastq(pipeconfig['path']['analysis']))
        for sl in symlink:
            os.symlink(*sl)
    except OSError:
        pass
    except:
        raise
    # add linked fastq
    inputfiles.append(ngsjob.linkedfastq())

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Start DRMAA session if available
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
try:
    from ruffus.drmaa_wrapper import run_job, error_drmaa_job
    import drmaa
except RuntimeError:  # could not load drmaa
    sge=None
else:
    sge = drmaa.Session()
    sge.initialize()

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Functions
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

def loadConfiguration(fi, configfile='.molpath'):
    # get path
    configpath = fi[:fi.rfind('/')]
    # load configuration from file
    try:
        fh = open('/'.join([configpath,configfile]))
    except:
        raise Exception("Could not load configuration (%s)" % configpath+configfile)
    else:
        ngsjob_config = NGSJob(json.load(fh))
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

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Pipeline tasks
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
from ruffus import *

@transform([inputfiles,trimmomatic],suffix('.fastq'),'_fastqc.html')
def fastQC(input_files,output_files):
    # load configurations
    p = loadConfiguration(input_files[0])
    # configure job
    job_name = 'fastqc'+p.RGSM()
    cmd = " ".join([ "fastqc", "--noextract", "--outdir="+p.wd(), ' '.join(input_files) ])
    # run on cluster
    if sge:
        try:
            stdout_res, stderr_res  = run_job(cmd,
                                                job_name          = job_name,
                                                logger            = logger,
                                                drmaa_session     = sge[0],
                                                run_locally       = False,
                                                job_other_options = "-q short.q",
                                                job_script_directory = "test_dir",
                                                job_environment={ 'BASH_ENV' : '~/.bashrc' },
                                                retain_job_scripts = True,
                                                working_directory = wd )
        except error_drmaa_job as err:
            raise Exception("\n".join(map(str,["Failed to run:",cmd,err,stdout_res,stderr_res])))
    # fallback
    else:
        run_cmd(cmd)

# this filters the fastq
@posttask (touch_file( 'trimming.completed' ))
@transform(inputfiles,
            formatter("(?P<PAIR>[^\.]+)\.fastq", "(?P<PAIR>[^\.]+)\.fastq"),
            ["{PAIR[0]}.filtered.fastq", "{PAIR[0]}.unpaired.fastq",
             "{PAIR[1]}.filtered.fastq", "{PAIR[1]}.unpaired.fastq"])
def trimmomatic(input_files, output_files):
    # load configurations
    p = loadConfiguration(input_files[0])
    # configure job
    job_name = 'trimmomatic'+p.RGSM()
    cmd = ' '.join([
        pipeconfig['software']['java'],
        '-XX:ParallelGCThreads='+str(pipeconfig['resources']['trimmomatic']['cpu']),
        '-Xmx'+str(pipeconfig['resources']['trimmomatic']['java_mem'])+'g',
        '-jar', pipeconfig['software']['trimmomatic'],
        p.librarytype(),
        '-phred64',
        '-trimlog', p.RGSM()+'.trimming.log',
        '-threads', str(pipeconfig['resources']['trimmomatic']['cpu']),
        ' '.join(input_files),
        ' '.join(output_files),
        ' '.join(ngsconfig['analysis'][p.ngsAnalysis]['config']['trimmomatic']) ])
    wd = p.wd()
    # run on cluster
    if sge:
        try:
            stdout_res, stderr_res  = run_job(cmd,
                                                job_name          = job_name,
                                                logger            = logger,
                                                drmaa_session     = sge[0],
                                                run_locally       = False,
                                                job_other_options = "-q short.q",
                                                job_script_directory = "test_dir",
                                                job_environment={ 'BASH_ENV' : '~/.bashrc' },
                                                retain_job_scripts = True,
                                                working_directory = wd )
        except error_drmaa_job as err:
            raise Exception("\n".join(map(str,["Failed to run:",cmd,err,stdout_res,stderr_res])))
    # fallback
    else:
        run_cmd(cmd)

# ALIGNER
### @collate(animals, regex(r"(.+)\.(.+)\.animal"),  r"\2.results")

@posttask (touch_file( 'alignment.completed' ))
@transform(inputfiles,
            formatter("(?P<PAIR>[^\.]+)\.fastq"),
            "{PAIR[0]}.filtered.fastq")
def alignment(input_files, output_files):
    print input_files
    print output_files
    sys.exit()

    # load configurations
    p = loadConfiguration(input_files[0])
    # configure job
    job_name = 'trimmomatic'+p.RGSM()
    cmd = ' '.join([
        pipeconfig['software']['java'],
        '-XX:ParallelGCThreads='+str(pipeconfig['resources']['trimmomatic']['cpu']),
        '-Xmx'+str(pipeconfig['resources']['trimmomatic']['java_mem'])+'g',
        '-jar', pipeconfig['software']['trimmomatic'],
        p.librarytype(),
        '-phred64',
        '-trimlog', p.RGSM()+'.trimming.log',
        '-threads', str(pipeconfig['resources']['trimmomatic']['cpu']),
        ' '.join(input_files),
        ' '.join(output_files),
        ' '.join(ngsconfig['analysis'][p.ngsAnalysis]['config']['trimmomatic']) ])
    wd = p.wd()
    # run on cluster
    if sge:
        try:
            stdout_res, stderr_res  = run_job(cmd,
                                                job_name          = job_name,
                                                logger            = logger,
                                                drmaa_session     = sge[0],
                                                run_locally       = False,
                                                job_other_options = "-q short.q",
                                                job_script_directory = "test_dir",
                                                job_environment={ 'BASH_ENV' : '~/.bashrc' },
                                                retain_job_scripts = True,
                                                working_directory = wd )
        except error_drmaa_job as err:
            raise Exception("\n".join(map(str,["Failed to run:",cmd,err,stdout_res,stderr_res])))
    # fallback
    else:
        run_cmd(cmd)


'''
@transform(trimmomatic,
    suffix(".filtered.fastq"),
    ".novoalign.sam")
def novoalign(input_files, output_files):
    print input_files, output_files
    return

@transform(trimmomatic,
    suffix(".filtered.fastq"),
    ".stampy.sam")
def stampy(input_files, output_files):
    print input_files, output_files
    return

# 3. sam2bam(p,sge):
# 4. sortSam(p, sge):
# 5. addReplaceReadGroups(p, sge):
# 6. MarkDuplicates(p, sge):
# 7. RealignerTargetCreator
# 8. IndelRealigner
# 9. BaseRecalibrator before recal
# 10. PrintReads_BQSR = QUAL SCORE RECALIBRATION
# 11. BaseRecalibrator after recal
# 12. AnalyzeCovariates before & after recal
# 13. HaplotypeCaller
# 14. VCFtoolsSiteFilter_on_HaplotypeCaller_Output
# 15. UnifiedGenotyper
# 16. VCFtoolsSiteFilter_on_UnifiedGenotyper_Output
# 17. BedTools_DepthOfCoverage
# 18. CollectMultipleMetrics
# 19. Table_Annovar
'''

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Print list of tasks
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
if options.just_print:
    pipeline_printout(sys.stdout, [trimmomatic], verbose=options.verbose)


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Print flowchart
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
elif options.flowchart:
    # use file extension for output format
    output_format = os.path.splitext(options.flowchart)[1][1:]
    pipeline_printout_graph (open(options.flowchart, "w"),
                             output_format,
                             [trimmomatic],
                             no_key_legend = True)
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Run Pipeline
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
else:
    pipeline_run([trimmomatic], forcedtorun_tasks = [trimmomatic], multiprocess = options.jobs, logger = logger, verbose=options.verbose)

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Cleanly end drmaa session
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
if sge:
    sge.exit()


