#!/usr/bin/env python

__doc__=="""
#############################################################################################
# NGS pipeline for molecular diagnostics (can be dockerized)                                #
# -- Organisation: KCL/SLaM/NHS/Viapath                                                     #
# -- Date: 19/08/2014                                                                       #
#############################################################################################
"""
__name__ = 'ngsEasy'
__author__ = "David Brawand"
__credits__ = ['Stephen Newhouse', 'Amos Folarin', 'Aditi Gulati']
__license__ = "LGPL"
__version__ = "0.9"
__maintainer__ = "David Brawand"
__email__ = "dbrawand@nhs.net, stephen.j.newhouse@gmail.com, amosfolarin@gmail.com, aditi.gulati@nhs.net"
__status__ = "Development"  # ["Prototype", "Development",  "Production"]


switch_WGS = True
switch_TAS = True

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   basic imports and paths
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
parser = OptionParser(version="%prog 1.0", usage=usage)

#   configuration files
parser.add_option("-p", dest="pipeconfig", default='pipeline_config.json',metavar="STRING", \
    help="Pipeline parameters (pipeline_config.json)")
parser.add_option("-n", dest="ngsconfig", default='ngs_config.json',metavar="STRING", \
    help="Workflow configuration (ngs_config.json)")

#   general options: verbosity / logging
parser.add_option("-v", dest="verbose", action="count", default=1, \
    help="verbosity (for ruffus messages)"
     "(eg. -vvv)")
parser.add_option("-l", dest="logfile", default=None, metavar="STRING", \
    help="logfile (default STDERR only)")
parser.add_option("-d" ,dest="debug", action='store_true', default=False,\
    help="display DEBUG messages and temporarily disable method checksum flow control")

#   pipeline options
parser.add_option("-c", dest="cleanup", default=False, action="store_true", \
    help="cleanup intermediary files (truncate to zero)")
parser.add_option("-t", dest="targettasks", default=None, metavar="task1,task2,...", \
    help="target tasks (run all tasks by default)")
parser.add_option("-f", dest="forcedtasks", default=None, metavar="task1,task2,...", \
    help="forced tasks (updates dependent tasks)")
parser.add_option("--touch", dest="touch", default=False, action="store_true", \
    help="just touch")
parser.add_option("--flowchart", dest="flowchart", metavar="FILE", type="string", \
                  help="Print flowchart of the pipeline to FILE. Flowchart format "
                       "depends on extension. Alternatives include ('.dot', '.jpg', "
                       "'*.svg', '*.png' etc). Formats other than '.dot' require "
                       "the dot program to be installed (http://www.graphviz.org/).")
parser.add_option("--just_print", dest="just_print", action="store_true", default=False, \
                    help="Only print a trace (description) of the pipeline. "
                         " The level of detail is set by --verbose.")
# Run options
parser.add_option("--jobs", dest="jobs", default=1, metavar="INT", type="int", \
    help="Specifies the number of jobs (operations) to run in parallel.")
parser.add_option("--local", dest="runlocal", default=False, action="store_true", \
    help="Force local execution (via DRMAA, if available)")

(options, args) = parser.parse_args()

if len(args) == 0:
    parser.print_help()
    parser.error("\n\tno sample file provided")

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Logger
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
import logging

# get logger and formatter
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
#formatter = logging.Formatter('%(asctime)s %(name)-10s %(levelname)-8s %(message)s')
formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s')

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

from ngsEasy_helpers import tail, zeroFile, parsePicard, flatten, timejob, githash

inputfiles = []
cfg_file = '.molpath'

# print ngsjob data in readable format
for i, ngsjob in enumerate(ngsjobs):
    targetDir = '/'.join([ pipeconfig["path"]["analysis"], ngsjob.wd() ])
    configpath = targetDir+'/'+cfg_file
    # build configfile and setup directories
    if not os.path.isfile(configpath):
        logger.info('Setting up '+str(ngsjob))

        # store GIT revision hash and pipeline configuration (subset)
        ngsjob.githash = githash()
        ngsjob._pipeconfig = {k: pipeconfig.get(k, None) for k in ('reference', 'software')}

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

        # get info from FastQ names (written to config file)
        for fq in ngsjob.fastq(pipeconfig['path']['fastq']):
            match_new = re.match(r".*?/([a-zA-Z0-9-.]+)_([^_/]+)(?:_[CAGTN]+)?_L([0-9]+)_R(1|2).fastq",fq)
            ## sample_name sameplenumber/run_id (barcode) lane pair ???
            if match_new:
                ngsjob._sample = match_new.group(1) if not ngsjob._sample else ','.join([ngsjob._sample, match_new.group(1)])
                ngsjob._run_id = match_new.group(2) if not ngsjob._run_id else ','.join([ngsjob._run_id, match_new.group(2)])
                ngsjob._lane = match_new.group(3) if not ngsjob._lane else ','.join([ngsjob._lane, match_new.group(3)])
                ngsjob._pair = match_new.group(4) if not ngsjob._pair else ','.join([ngsjob._pair, match_new.group(4)])
            else:
                logger.warn("Unable to parse name of fastq file %s" % fq)

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
    logger.warn("no DRMAA found -> runs locally")
    options.multiprocess = 1
else:
    if pipeconfig['docker']['active']:
        sge = None
        sys.exit("### DOCKERIZED DISPATCH NOT IMPLEMENTED YET ###")
    elif pipeconfig['switch']['forcelocal']:
        sge = None
        logger.info("Forced local run (this may take significantly longer)")
    else:
        sge = drmaa.Session()
        sge.initialize()
        logger.info("using DRMAA for job dispatch")

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Functions
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

'''load configuration file in directory'''
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

'''controlled local execution of command'''
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
'''run task on a computing grid (using ruffus.drmaa wrapper)'''
def run_sge(cmd_str, jobname, fullwd, cpu=1, mem=2, runlocally=False, scriptdir="tmp" , sge_env={ 'BASH_ENV' : '~/.bashrc' }):
    # get options
    opt = []
    opt.append('-o '+fullwd+'/sge.out')
    opt.append('-e '+fullwd+'/sge.out')
    for k,v in pipeconfig['resources']['gridengine'].items():
        opt.append(' '.join([k,v]))
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
        raise Exception("\n".join(map(str,["Failed to run:",cmd_str,'with error',err])))

'''writes shell script to execute in a docker container'''
def dockerDispatch(cmd_str, run_container):
    logger.debug("Running \"%s\" in container \"%s\"" % (cmd_str, run_container))
    return


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   get input files
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
for i, ngsjob in enumerate(ngsjobs):
    logger.info("Reading configuration from %s" % targetDir)
    targetDir = '/'.join([ pipeconfig["path"]["analysis"], ngsjob.wd() ])
    inputfiles.append(loadConfiguration(targetDir).inputfiles())

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   flowchart styles
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
style_start = {'shape': "box"}
style_collect = {  'shape': "hexagon" }
style_normal = { 'shape': "box" }
style_end = {'shape': "box"}

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Pipeline tasks
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
from ruffus import *
import inspect  # to get the functions name inside
from collections import Counter, defaultdict
from numpy import median

#--------------------
# just count reads
#--------------------
@graphviz(label_prefix="FastQ check\n", **style_start)
@transform(inputfiles, suffix("_1.fastq"), '.check')
@timejob(logger)
def ngsEasy_checkPaired(input_files,output_file):
    '''
    checks if first and last read are paired, assumes miSeq format
    '''
    # get mates
    with open(input_files[0], 'r') as f:
        firstmate =  (f.readline()[1:].split(' ')[0], tail(f)[1:].split('\n')[0].split(' ')[0])
    with open(input_files[1], 'r') as f:
        secondmate =  (f.readline()[1:].split(' ')[0], tail(f)[1:].split('\n')[0].split(' ')[0])
    #print firstmate
    #print secondmate
    # check if equal
    if firstmate[0] == secondmate[0] and firstmate[1] == secondmate[1]:
        # write info to file
        first = re.match(r"([^:]+):(\d+):([^:]+):(\d+):(\d+):\d+:\d+" ,firstmate[0])
        last = re.match(r"([^:]+):(\d+):([^:]+):(\d+):(\d+):\d+:\d+" ,firstmate[1])
        with open(output_file,'w') as ofh:
            try:
                ofh.write('InstrumentId\t'+first.group(1)+'\n')
                ofh.write('runId\t'+first.group(2)+'\n')
                ofh.write('flowcellId\t'+first.group(3)+'\n')
                ofh.write('flowcellLane'+'\t'+first.group(4)+'\n')
                ofh.write('firstTile'+'\t'+first.group(5)+'\n')
                ofh.write('lastTile'+'\t'+last.group(5)+'\n')
            except:
                raise
    else:
        logger.debug('firstMate  %s %s' % firstmate)
        logger.debug('secondMate %s %s' % secondmate)
        logger.critical('FastQ are not paired or incomplete. Please check.')
        sys.exit(1)

#--------------------
# FASTQC prefiltering
#--------------------
@graphviz(label_prefix="inital FastQC\n", **style_normal)
@follows(ngsEasy_checkPaired)
@transform(inputfiles, formatter("(?P<PAIR>[^\.]+)\.fastq", "(?P<PAIR>[^\.]+)\.fastq"),
            ["{PAIR[0]}_fastqc.zip", "{PAIR[1]}_fastqc.zip"])
@timejob(logger)
def ngsEasy_prefilterFastQC(input_files,output_files):
    # load configurations and name job
    p = loadConfiguration(input_files[0][:input_files[0].rfind('/')])
    cmd = " ".join([ pipeconfig['software']['fastqc'], "--noextract"] + input_files)
    # # run on cluster
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=1, mem=2,runlocally=options.runlocal)
    else:
        run_cmd(cmd)

#--------------------
# Adapter trimming
#--------------------
#@posttask (touch_file( 'trimming.completed' ))
#@jobs_limit(8,'iolimit')  # further limit jobs with high I/O
@graphviz(label_prefix="Adapter trimming\n", **style_normal)
@follows(ngsEasy_checkPaired)
@transform(inputfiles,
            formatter("(?P<PAIR>[^\.]+)\.fastq", "(?P<PAIR>[^\.]+)\.fastq"),
            ["{PAIR[0]}.filtered.fastq", "{PAIR[1]}.filtered.fastq"])
@timejob(logger)
def ngsEasy_trimmomatic(input_files, output_files):
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
        '-trimlog', input_files[0].replace('_1.fastq','.trimlog'),
        '-threads', str(pipeconfig['resources']['trimmomatic']['cpu']),
        ' '.join(input_files),
        ' '.join(filtered),
        ' '.join(ngsconfig['ngsanalysis'][p.ngsAnalysis]['trimmomatic']) ])
    # run on cluster
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=pipeconfig['resources']['trimmomatic']['cpu'], mem=pipeconfig['resources']['trimmomatic']['mem'] , runlocally=options.runlocal)
    # fallback
    else:
        run_cmd(cmd)
    # cleanup (discards unpaired reads)
    if options.cleanup:
        for filteredfile in filtered:
            if filteredfile not in output_files:
                zeroFile(filteredfile)

#--------------------
# FASTQC postfilter
#--------------------
@graphviz(label_prefix="Final FastQC\n", **style_normal)
@transform(ngsEasy_trimmomatic, formatter("(?P<PAIR>[^\.]+\.filtered)\.fastq", "(?P<PAIR>[^\.]+\.filtered)\.fastq"),
            ["{PAIR[0]}_fastqc.zip", "{PAIR[1]}_fastqc.zip"])
@timejob(logger)
def ngsEasy_postfilterFastQC(input_files,output_files):
    # load configurations
    p = loadConfiguration(input_files[0][:input_files[0].rfind('/')])
    cmd = " ".join([ pipeconfig['software']['fastqc'], "--noextract"] + input_files)
    # run on cluster
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=1, mem=2,runlocally=options.runlocal)
    # fallback
    else:
        run_cmd(cmd)

#--------------------
# read QC report
#--------------------
###### ADD MERGE/COLLATE OF FASTQC and gerneate QC report
@graphviz(label_prefix="Read QC\n", **style_collect)
@follows(ngsEasy_postfilterFastQC)
@collate(ngsEasy_prefilterFastQC, regex(r'(.+)_1_fastqc\.zip'), r'\1.readQC')
@timejob(logger)
def ngsEasy_readQC(input_files,output_file):
    p = loadConfiguration(input_files[0][0][:input_files[0][0].rfind('/')])
    # check trimlog
    try:
        filename =  input_files[0][0]
        m = re.match('/(.+\/)*([^_]+)',filename)
        assert m
    except:
        logger.warn('Did not find trimming log (pass)')
        with open(output_file,'a'): pass
    else:
        re_trimlog = re.compile('\S+\s+(\d)\S+\s+(\d+)\s+\d+\s+\d+\s+\d+$')
        success = Counter()
        lengths = defaultdict(list)
        trimlog = '/'+m.group(1)+m.group(2)+'.trimlog'
        for l in open(trimlog):
            m = re_trimlog.match(l)
            if int(m.group(2)) != 0:
                success[m.group(1)] += 1
            lengths[m.group(1)].append(int(m.group(2)))
        # check how many reads remain (check median read length)
        with open(output_file,'w') as rqc:
            for k in lengths.keys():
                _median = median(lengths[k])
                print >> rqc, k, 'median:'+str(_median),
                if _median >= ngsconfig['ngsanalysis'][p.ngsAnalysis]['readQC']['medianlength']:
                    print >> rqc, "PASS"
                else:
                    print >> rqc, "FAIL"
                print >> rqc, k, 'count:'+str(success[k]),
                if success[k] >= ngsconfig['ngsanalysis'][p.ngsAnalysis]['readQC']['count']:
                    print >> rqc, "PASS"
                else:
                    print >> rqc, "FAIL"
                # abort trap
                try:
                    assert success[k] > 0
                except AssertionError:
                    logger.error("ABORTING: NO READS REMAIN AFTER FILTERING!")
                    sys.exit(1)
                except:
                    raise

#--------------------
# ALIGNMENT
#--------------------
@follows(ngsEasy_readQC)
@graphviz(label_prefix="Align reads\n", **style_normal)
@transform(ngsEasy_trimmomatic, formatter("(?P<PAIR>[^\.]+)_1\.filtered\.fastq"), "{PAIR[0]}.sam")
@timejob(logger)
def ngsEasy_alignment(input_files, output_file):
    p = loadConfiguration(input_files[0][:input_files[0].rfind('/')])
    # check aligner
    aligner = ngsconfig['ngsanalysis'][p.ngsAnalysis]['aligner']
    try:
        assert aligner in ['bwa','bowtie2','stampy','bowtie','bowtie2','novoalign']
    except:
        logger.warn("aligner %s not available. Falling back to bwa." % aligner)
        aligner = 'bwa'
    # GE resources
    try:
        ge_cpu = pipeconfig['resources'][aligner]['cpu']
        ge_mem = pipeconfig['resources'][aligner]['mem']
    except:
        logger.warn("Aligner (%s) cpu/mem set to defaults" % aligner)
        ge_cpu = 4
        ge_mem = 8
    # configure job
    if aligner.startswith('bwa'):
        cmd = " ".join([
            pipeconfig['software']['bwa'],
            'mem', '-M',
            '-t '+str(pipeconfig['resources']['bwa']['cpu']),
            pipeconfig['reference']['genome']['bwaindex'],
            ' '.join(input_files),
            '>', output_file ])
    elif aligner.startswith('stampy'):
        bwa = output_file.replace('.sam','.bwa.sam')
        cmd = " ".join([
            pipeconfig['software']['bwa'],
            'mem', '-M',
            '-t '+ str(pipeconfig['resources']['bwa']['cpu']),
            pipeconfig['reference']['genome']['sequence'],
            ' '.join(input_files),
            #'|', pipeconfig['software']['bwa'], 'view -Sb -',
            '>', bwa,
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
    elif aligner.startswith('novo'):
        stat_file = output_file.replace('.sam','.stat')
        cmd = " ".join([
            pipeconfig['software']['novoalign'],
            '-d', pipeconfig['reference']['genome']['novoindex'],
            '-f', ' '.join(input_files),
            '-F STDFQ --Q2Off --3Prime -g 40 -x 6 -r All -i PE 300,150',
            '-c', pipeconfig['resources']['novoalign']['cpu'],
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
            cpu=ge_cpu, mem=ge_mem, runlocally=options.runlocal)
    else:
        run_cmd(cmd)
    # cleanup of intermediary files (zero files)
    if options.cleanup:
        for infile in input_files:
            zeroFile(infile)

#--------------------
# process SAM -> BAM
#--------------------
@transform(ngsEasy_alignment, suffix('.sam'), '.bam')
@graphviz(label_prefix="SAM -> BAM\n", **style_normal)
@timejob(logger)
def ngsEasy_sam2bam(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    cmd = " ".join([ pipeconfig['software']['samtools'], "view", "-Shb", '-o', output_file, input_file ])
    # run on cluster
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=1, mem=2,runlocally=options.runlocal)
    else:
        run_cmd(cmd)
    # cleanup
    if options.cleanup:
        zeroFile(input_file)


@transform(ngsEasy_sam2bam, suffix('.bam'), '.sorted.bam')
@graphviz(label_prefix="Sort aligned reads\n", **style_normal)
@timejob(logger)
def ngsEasy_sortSam(input_file,output_file):
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
            cpu=1, mem=2,runlocally=options.runlocal)
    else:
        run_cmd(cmd)
    # cleanup
    if options.cleanup:
        zeroFile(input_file)


@transform(ngsEasy_sortSam, suffix('.bam'), '.bai')
@graphviz(label_prefix="Index reads\n", **style_normal)
@timejob(logger)
def ngsEasy_indexBam(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    cmd = " ".join([ pipeconfig['software']['samtools'], "index", input_file, output_file ])
    # run on cluster
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=1, mem=2,runlocally=options.runlocal)
    else:
        run_cmd(cmd)

#--------------------
# add/replace ReadGroup
#--------------------
@follows(ngsEasy_indexBam)
@transform(ngsEasy_sortSam, suffix('.sorted.bam'), '.addrg.bam')
@graphviz(label_prefix="set Read Groups\n", **style_normal)
@timejob(logger)
def ngsEasy_addReplaceReadGroups(input_file,output_file):
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
            cpu=1, mem=pipeconfig['resources']['picard']['mem'],runlocally=options.runlocal)
    else:
        run_cmd(cmd)
    # cleanup
    if options.cleanup:
        zeroFile(input_file)

#--------------------
# mark/remove Duplicates
#--------------------
@transform(ngsEasy_addReplaceReadGroups, suffix('.addrg.bam'), '.dupemk.bam')
@graphviz(label_prefix="Mark PCR duplicates\n", **style_normal)
@timejob(logger)
def ngsEasy_markDuplicates(input_file,output_file):
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
            cpu=1, mem=2,runlocally=options.runlocal)
    else:
        run_cmd(cmd)
    # cleanup
    if options.cleanup:
        zeroFile(input_file)

#--------------------
# Alignment QC
#--------------------
@graphviz(label_prefix="Collect Alignment metrics\n", **style_normal)
@transform(ngsEasy_markDuplicates, suffix('.dupemk.bam'), \
    ['.metrics.alignment_summary_metrics',
    '.metrics.insert_size_metrics',
    '.metrics.quality_distribution_metrics',
    '.metrics.quality_distribution.pdf',
    '.metrics.quality_by_cycle.pdf'])
@timejob(logger)
def ngsEasy_collectMultipleMetrics(input_file,output_files):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    cmd = " ".join([
        pipeconfig['software']['java'],
        '-XX:ParallelGCThreads='+str(pipeconfig['resources']['picard']['cpu']),
        '-Xmx'+str(pipeconfig['resources']['picard']['java_mem'])+'g',
        '-jar', pipeconfig['software']['picard'].rstrip('/')+'/CollectMultipleMetrics.jar',
        'TMP_DIR='+'/'.join([pipeconfig['path']['analysis'], p.wd(), 'tmp']),
        'VALIDATION_STRINGENCY=SILENT',
        'MAX_RECORDS_IN_RAM=100000',
        'INPUT='+input_file,
        'OUTPUT='+output_files[0][:output_files[0].rfind('.')],
        'REFERENCE_SEQUENCE='+pipeconfig['reference']['genome']['sequence'],
        #NOT NEEDED (in inital metrics)'PROGRAM=CollectAlignmentSummaryMetrics',
        'PROGRAM=CollectInsertSizeMetrics',
        'PROGRAM=QualityScoreDistribution',
        'PROGRAM=MeanQualityByCycle'
        ])
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=1, mem=2,runlocally=options.runlocal)
    else:
        run_cmd(cmd)

@active_if(switch_TAS)
@graphviz(label_prefix="Calculate target metrics\n", **style_normal)
@follows(ngsEasy_collectMultipleMetrics)
@transform(ngsEasy_markDuplicates, suffix('.dupemk.bam'), ['.metrics.TAS','.metrics.TAS.coverage'])
@timejob(logger)
def ngsEasy_collectTargetedPcrMetrics(input_file,output_files):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    ampliconIntervals = ngsconfig['ngsanalysis'][p.ngsAnalysis]['tas']['bait_intervals']
    targetIntervals = ngsconfig['ngsanalysis'][p.ngsAnalysis]['tas']['target_intervals']
    # command
    cmd = " ".join([
        pipeconfig['software']['java'],
        '-XX:ParallelGCThreads='+str(pipeconfig['resources']['picard']['cpu']),
        '-Xmx'+str(pipeconfig['resources']['picard']['java_mem'])+'g',
        '-jar', pipeconfig['software']['picard'].rstrip('/')+'/CollectTargetedPcrMetrics.jar',
        'TMP_DIR='+'/'.join([pipeconfig['path']['analysis'], p.wd(), 'tmp']),
        'VALIDATION_STRINGENCY=SILENT',
        'MAX_RECORDS_IN_RAM=100000',
        'REFERENCE_SEQUENCE='+pipeconfig['reference']['genome']['sequence'],
        'AMPLICON_INTERVALS='+ampliconIntervals,
        '{targets}'.format(targets='TARGET_INTERVALS='+targetIntervals if targetIntervals else ""),
        'PER_TARGET_COVERAGE='+output_files[1],
        'METRIC_ACCUMULATION_LEVEL='+'ALL_READS', # sample, library, read_group
        'CUSTOM_AMPLICON_SET_NAME='+p.ngsType,
        'INPUT='+input_file,
        'OUTPUT='+output_files[0]
        ])
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=1, mem=2,runlocally=options.runlocal)
    else:
        run_cmd(cmd)

@active_if(switch_WGS)
@graphviz(label_prefix="Calculate WGS metrics\n", **style_normal)
@transform(ngsEasy_markDuplicates, suffix('.dupemk.bam'), '.metrics.WGS')
@timejob(logger)
def ngsEasy_collectWgsMetrics(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    cmd = " ".join([
        pipeconfig['software']['java'],
        '-XX:ParallelGCThreads='+str(pipeconfig['resources']['picard']['cpu']),
        '-Xmx'+str(pipeconfig['resources']['picard']['java_mem'])+'g',
        '-jar', pipeconfig['software']['picard'].rstrip('/')+'/CollectWgsMetrics.jar',
        'TMP_DIR='+'/'.join([pipeconfig['path']['analysis'], p.wd(), 'tmp']),
        'VALIDATION_STRINGENCY=SILENT',
        'MAX_RECORDS_IN_RAM=100000',
        'INPUT='+input_file,
        'OUTPUT='+output_file,
        'REFERENCE_SEQUENCE='+pipeconfig['reference']['genome']['sequence'],
        'MINIMUM_MAPPING_QUALITY='+str(20),
        'MINIMUM_BASE_QUALITY='+str(20),
        'COVERAGE_CAP='+str(1000)
        ])
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=1, mem=2,runlocally=options.runlocal)
    else:
        run_cmd(cmd)

#--------------------
# alignment QC report
#--------------------
'''writes picard results into a single flat file'''
@timejob(logger)
@graphviz(label_prefix="Alignment QC\n", **style_collect)
@collate([ngsEasy_collectMultipleMetrics, ngsEasy_collectTargetedPcrMetrics, ngsEasy_collectWgsMetrics],
    regex(r'(.+)\.metrics\.[^\.]+'),
    r'\1.summarymetrics')
@timejob(logger)
def ngsEasy_alignmentQC(input_files,output_file):
    with open(output_file,'w') as outfh:
        for infile in flatten(input_files):
            if infile.endswith('pdf') or infile.endswith('coverage'):
                continue
            with open(infile) as infh:
                fields, ordering, postfix = parsePicard(infh)
                # write summary to file
                for n in sorted(ordering.keys()):
                    k = ordering[n]
                    v = fields[k]
                    if len(v) == 0:
                        continue
                    for i,p in enumerate(postfix):
                        try:
                            float(v[i])
                        except ValueError:
                            pass
                        else:
                            print >> outfh, '{:<30} {:<16} {:<26} {}'.format(infile[infile.rfind('.')+1:], postfix[i], k, v[i])

'''
__FUTURE__DEVELOPMENT__
#upload to QC database
def storeQC(input_file, output_file):
    # save in sql database
    pass
'''

#--------------------
# GATK indelRealign
#--------------------
@graphviz(label_prefix="Create realignment targets\n", **style_normal)
@transform(ngsEasy_markDuplicates, suffix('.dupemk.bam'), '.dupemk.intervals')
@timejob(logger)
def ngsEasy_realignerTargetCreator(input_file,output_file):
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
            cpu=pipeconfig['resources']['gatk']['RTC']['cpu'], mem=pipeconfig['resources']['gatk']['RTC']['mem'],runlocally=options.runlocal)
    else:
        run_cmd(cmd)

@graphviz(label_prefix="Realign indels\n", **style_normal)
@collate([ngsEasy_markDuplicates,ngsEasy_realignerTargetCreator], regex(r'(.+)\.dupemk\.[^\.]+'), r'\1.realn.bam')
@timejob(logger)
def ngsEasy_indelRealigner(input_files,output_file):
    p = loadConfiguration(input_files[0][:input_files[0].rfind('/')])
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
            cpu=pipeconfig['resources']['gatk']['IR']['cpu'], mem=ipeconfig['resources']['gatk']['IR']['mem'],runlocally=options.runlocal)
    else:
        run_cmd(cmd)
    # cleanup
    if options.cleanup:
        for infile in input_files:
            zeroFile(infile)

@graphviz(label_prefix="Index aligned reads\n", **style_normal)
@transform(ngsEasy_indelRealigner, suffix('.bam'), '.bai')
@timejob(logger)
def ngsEasy_indexRealignedBam(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    cmd = " ".join([ pipeconfig['software']['samtools'], "index", input_file, output_file ])
    # run on cluster
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=1, mem=2,runlocally=options.runlocal)
    else:
        run_cmd(cmd)

#--------------------
# GATK BaseRecalib
#--------------------

@follows(ngsEasy_indexRealignedBam)
@graphviz(label_prefix="Calculate base calibration\n", **style_normal)
@transform(ngsEasy_indelRealigner, suffix('.realn.bam'), '.realn.recalibrationtable')
@timejob(logger)
def ngsEasy_firstBaseRecalibrator(input_file,output_file):
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
            cpu=pipeconfig['resources']['gatk']['BR']['cpu'], mem=pipeconfig['resources']['gatk']['BR']['mem'],runlocally=options.runlocal)
    else:
        run_cmd(cmd)

@follows(ngsEasy_indexRealignedBam)
@graphviz(label_prefix="Recalibrate aligned reads\n", **style_normal)
@collate([ngsEasy_indelRealigner,ngsEasy_firstBaseRecalibrator], regex(r'(.+)\.realn\.[^\.]+'), r'\1.recal.bam')
@timejob(logger)
def ngsEasy_recalibrateBam(input_files,output_file):
    p = loadConfiguration(input_files[0][:input_files[0].rfind('/')])
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
        '-I', input_files[0],
        '-o', output_file
        ])
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=pipeconfig['resources']['gatk']['PR']['cpu'], mem=pipeconfig['resources']['gatk']['PR']['mem'],runlocally=options.runlocal)
    else:
        run_cmd(cmd)
    # cleanup
    if options.cleanup:
        for infile in input_files:
            if not infile.endswith('recalibrationtable'):
                zeroFile(infile)

@transform(ngsEasy_recalibrateBam, suffix('.bam'), '.bai')
@graphviz(label_prefix="Index recalibrated reads\n", **style_normal)
@timejob(logger)
def ngsEasy_indexRecalibratedBam(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    cmd = " ".join([ pipeconfig['software']['samtools'], "index", input_file, output_file ])
    # run on cluster
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=1, mem=2,runlocally=options.runlocal)
    else:
        run_cmd(cmd)

@follows(ngsEasy_indexRecalibratedBam)
@graphviz(label_prefix="Calculate base calibration\n", **style_normal)
@transform(ngsEasy_recalibrateBam, suffix('.recal.bam'), '.recal.recalibrationtable')
@timejob(logger)
def ngsEasy_finalBaseRecalibrator(input_file,output_file):
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
            cpu=pipeconfig['resources']['gatk']['BR']['cpu'], mem=pipeconfig['resources']['gatk']['BR']['mem'],runlocally=options.runlocal)
    else:
        run_cmd(cmd)

#@collate([ngsEasy_firstBaseRecalibrator,ngsEasy_finalBaseRecalibrator], regex(r'(.+)\.(re...)\.recalibrationtable'), [r'\1.covariates.csv',r'\1.covariates.pdf'])
@graphviz(label_prefix="Calculate recalibration covariates\n", **style_collect)
@collate([ngsEasy_firstBaseRecalibrator,ngsEasy_finalBaseRecalibrator], regex(r'(.+)\.(re...)\.recalibrationtable'), r'\1.covariates.csv')
@timejob(logger)
def ngsEasy_analyzeCovariates(input_files,output_file):
    p = loadConfiguration(input_files[0][:input_files[0].rfind('/')])
    cmd = " ".join([
        pipeconfig['software']['java'],
        '-Xmx'+str(pipeconfig['resources']['gatk']['AC']['java_mem'])+'g',
        '-Djava.io.tmpdir='+'/'.join([pipeconfig['path']['analysis'], p.wd(), 'tmp']),
        '-jar', pipeconfig['software']['gatk'],
        '-T', 'AnalyzeCovariates',
        '-R', pipeconfig['reference']['genome']['sequence'],
        '-before', input_files[0],
        '-after', input_files[1],
        #'-plots', output_files[1],
        '-csv', output_file
        ])
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=pipeconfig['resources']['gatk']['AC']['cpu'], mem=pipeconfig['resources']['gatk']['AC']['mem'],runlocally=options.runlocal)
    else:
        run_cmd(cmd)
    # cleanup
    if options.cleanup:
        for infile in input_files:
            zeroFile(infile)

#--------------------
# Haplotype Calling UG
#--------------------
@active_if(pipeconfig['switch']['UG'])
@follows(ngsEasy_indexRecalibratedBam)
@graphviz(label_prefix="Unified Genotyper\n", **style_normal)
@transform(ngsEasy_recalibrateBam, suffix('.recal.bam'), '.UG.vcf')
@timejob(logger)
def ngsEasy_unifiedGenotyper(input_file,output_file):
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
            cpu=pipeconfig['resources']['gatk']['UG']['cpu'], mem=pipeconfig['resources']['gatk']['UG']['mem'],runlocally=options.runlocal)
    else:
        run_cmd(cmd)

#--------------------
# Haplotype Calling HC
#--------------------
@active_if(pipeconfig['switch']['HC'])
@follows(ngsEasy_indexRecalibratedBam)
@graphviz(label_prefix="Haplotype Caller\n", **style_normal)
@transform(ngsEasy_recalibrateBam, suffix('.recal.bam'), '.HC.vcf')
@timejob(logger)
def ngsEasy_haplotypeCaller(input_file,output_file):
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
            cpu=pipeconfig['resources']['gatk']['HC']['cpu'], mem=pipeconfig['resources']['gatk']['HC']['mem'],runlocally=options.runlocal)
    else:
        run_cmd(cmd)

#--------------------
# mpileup variant calling
#--------------------
@active_if(pipeconfig['switch']['ST'])
@follows(ngsEasy_indexRecalibratedBam)
@graphviz(label_prefix="MPILEUP variant calling\n", **style_normal)
@transform(ngsEasy_recalibrateBam, suffix('.recal.bam'), '.ST.vcf')
@timejob(logger)
def ngsEasy_samtoolsMpileup(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    try:
        pipeconfig['software']['samtools1']
    except KeyError:
        logger.warn("Didn't find Samtools 1.0+, using legacy samtools")
        cmd = " ".join([
            pipeconfig['software']['samtools'],
            'mpileup',
            '-g',
            '-f', pipeconfig['reference']['genome']['sequence'],
            input_file,
            '|',
            pipeconfig['software']['bcftools'],
            'call',
            '-m',
            #'-v', # variants only
            '-O v',
            '-o', output_file
            ])
    else:
        cmd = " ".join([
            pipeconfig['software']['samtools1'],
            'mpileup',
            '-vu',  # for samtools-1.0 (which doesnt index bams currently)
            '-f', pipeconfig['reference']['genome']['sequence'],
            input_file,
            '|',
            pipeconfig['software']['bcftools'],
            'call',
            '-m',
            #'-v', # variants only
            '-O v',
            '-o', output_file
            ])
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=pipeconfig['resources']['samtools']['cpu'], mem=pipeconfig['resources']['samtools']['mem'],runlocally=options.runlocal)
    else:
        run_cmd(cmd)

#--------------------
# platypus variant calling
#--------------------
@active_if(pipeconfig['switch']['PL'])
@graphviz(label_prefix="Platypus variant calling\n", **style_normal)
@transform(ngsEasy_markDuplicates, suffix('.dupemk.bam'), '.PL.vcf')
@timejob(logger)
def ngsEasy_platypus(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    cmd = " ".join([
        'python',
        pipeconfig['software']['platypus'],
        'callVariants',
        '--bamFiles='+input_file,
        '--refFile='+pipeconfig['reference']['genome']['sequence'],
        #NAIVE CALLING'--source='+pipeconfig['reference']['snps'][0]+'.gz',
        #'--minPosterior='+str(0),
        '--outputRefCalls=1',
        '--output='+output_file,
        '--nCPU='+str(pipeconfig['resources']['platypus']['cpu']),
        '--minReads='+str(ngsconfig['ngsanalysis'][p.ngsAnalysis]['platypus']['minReads']),
        '--filterDuplicates='+str(ngsconfig['ngsanalysis'][p.ngsAnalysis]['platypus']['filterDuplicates']),
        '--trimOverlapping=1'
        ])
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=pipeconfig['resources']['platypus']['cpu'], mem=pipeconfig['resources']['platypus']['mem'],runlocally=options.runlocal)
    else:
        run_cmd(cmd)

#--------------------
# freebayes variant calling
#--------------------
@active_if(pipeconfig['switch']['FB'])
@graphviz(label_prefix="Freebayes variant calling\n", **style_normal)
@transform(ngsEasy_markDuplicates, suffix('.dupemk.bam'), '.FB.vcf')
@timejob(logger)
def ngsEasy_freebayes(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    cmd = " ".join([
        pipeconfig['software']['freebayes'],
        '-f', pipeconfig['reference']['genome']['sequence'],
        '--report-monomorphic',
        '-C', str(ngsconfig['ngsanalysis'][p.ngsAnalysis]['freebayes']['minObs']),
        '{usedup}'.format(usedup="--use-duplicate-reads" if ngsconfig['ngsanalysis'][p.ngsAnalysis]['freebayes']['useDuplicates'] else ""),
        '-v', output_file,
        input_file
        ])
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=pipeconfig['resources']['freebayes']['cpu'], mem=pipeconfig['resources']['freebayes']['mem'],runlocally=options.runlocal)
    else:
        run_cmd(cmd)

#--------------------
# compress and index VCF
#--------------------
@graphviz(label_prefix="Compress and Index variants\n", **style_normal)
@transform([ngsEasy_unifiedGenotyper,ngsEasy_haplotypeCaller,ngsEasy_samtoolsMpileup,ngsEasy_platypus,ngsEasy_freebayes],
    regex(r'(.+)\.vcf'), r'\1.vcf.gz.tbi')
@timejob(logger)
def ngsEasy_compressIndexVCF(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    # quality filtering
    cmds = [
        [ pipeconfig['software']['bgzip'], '-c', '-f', input_file, '>', input_file+'.gz' ],
        [ pipeconfig['software']['tabix'], input_file+'.gz' ]
            ]
    # run job
    for i, cmd in enumerate(cmds):
        if sge:
            run_sge(' '.join(cmd),
                jobname="_".join([inspect.stack()[0][3], str(i+1), p.RG('SM')]),
                fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
                cpu=1, mem=2,runlocally=options.runlocal)
        else:
            run_cmd(' '.join(cmd))

#--------------------
# add missing dbSNP IDs with vcf-tools
#--------------------
@graphviz(label_prefix="Adding dbSNP IDs\n", **style_normal)
@transform(ngsEasy_compressIndexVCF, regex(r'(.+)(\.vcf.gz).tbi'), r'\1.dbsnp\2')
@timejob(logger)
def ngsEasy_identifyVariants(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    # convert to bcf
    cmd = [
        pipeconfig['software']['bcftools'],
        'annotate',
        '-a', pipeconfig['reference']['snps'][0]+'.gz',
        '-c', 'CHROM,POS,REF,ALT,ID,-,-,-',
        '-O z',
        '-o', output_file,
        input_file.replace('.tbi','') ]
    # run job
    if sge:
        run_sge(' '.join(cmd),
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=1, mem=4,runlocally=options.runlocal)
    else:
        run_cmd(' '.join(cmd))

@graphviz(label_prefix="Index VCF\n", **style_normal)
@transform(ngsEasy_identifyVariants, suffix('.vcf.gz'), '.vcf.gz.tbi')
@timejob(logger)
def ngsEasy_indexVCF(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    # quality filtering
    cmd = [ pipeconfig['software']['tabix'], input_file ]
    # run job
    if sge:
        run_sge(' '.join(cmd),
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=1, mem=2,runlocally=options.runlocal)
    else:
        run_cmd(' '.join(cmd))

#--------------------
# generate vcf stats
#--------------------
@follows(ngsEasy_indexVCF)
@graphviz(label_prefix="Variant Statistics\n", **style_normal)
@transform(ngsEasy_identifyVariants, suffix('.vcf.gz'), '.stats')
@timejob(logger)
def ngsEasy_varStats(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    # filtering
    cmd = [
        pipeconfig['software']['bcftools'],
        'stats',
        '-F', pipeconfig['reference']['genome']['sequence'],
        '-s', '-',  # all samples
        input_file,
        '>', output_file
        ]
    # run job
    if sge:
        run_sge(' '.join(cmd),
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(), cpu=1, mem=2,runlocally=options.runlocal)
    else:
        run_cmd(' '.join(cmd))

#--------------------
# filter all VCFs (quality)
#--------------------
@follows(ngsEasy_indexVCF)
@graphviz(label_prefix="Filter variants\n", **style_normal)
@transform(ngsEasy_identifyVariants, suffix('.dbsnp.vcf.gz'), '.filtered.vcf.gz')
@timejob(logger)
def ngsEasy_siteFilter(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    # filter
    filterstring = []
    #filterstring += ['-s', 'LOWQUAL']  # softfilter
    m = re.match(r'.+\.(..)\.dbsnp\.vcf\.gz', input_file)  # for caller dependent filtering
    try:
        assert m
    except:
        logger.error('Cannot guess variant caller')
        raise
    else:
        varCaller = m.group(1)
    # get exclusion (filtering) string
    try:
        assert len(ngsconfig['ngsanalysis'][p.ngsAnalysis]['variantfilter'][varCaller])>0
    except AssertionError:
        pass  # empty exclusion string
        logger.warn('empty filtering string for varCaller %s' % varCaller)
    except IndexError:
        logger.error('no filter string available for %s' % varCaller)  # no filter string for caller (configuration error)
        sys.exit(1)
    else:
        filterstring += ['-e', r"'"+ngsconfig['ngsanalysis'][p.ngsAnalysis]['variantfilter'][varCaller]+r"'"]
    # filtering
    cmd = [
        pipeconfig['software']['bcftools'],
        'filter',
        '-O', 'z',
        '-o', output_file,
        ' '.join(filterstring),
        input_file ]
    # run job
    if sge:
        run_sge(' '.join(cmd),
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=1, mem=2,runlocally=options.runlocal)
    else:
        run_cmd(' '.join(cmd))

#--------------------
# assess total number of genotyped dbSNP ids
#--------------------


#--------------------
# reduce VCF to targets and variant sites
#--------------------
@graphviz(label_prefix="Filter Regions\n", **style_normal)
@transform(ngsEasy_siteFilter, suffix('.filtered.vcf.gz'), '.targets.vcf')
@timejob(logger)
def ngsEasy_regionFilter(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    # reduce to target regions of targeted sequencing or exome
    if ngsconfig['ngstype'][p.ngsType]['RGLB'].upper() in ['TAS','EXO']:
        cmd = [
            pipeconfig['software']['bedtools'],
            'intersect',
            '-header',
            '-a', input_file,
            '-b', ngsconfig['ngsanalysis'][p.ngsAnalysis]['tas']['covered_bed'],
            '>', output_file
            ]
    else:
        logger.info('No region filtering for sample %s' % p.RG('SM'))
        cmd = [
            'cat', input_file,
            '>', output_file
        ]
    # run job
    if sge:
        run_sge(' '.join(cmd),
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=1, mem=2,runlocally=options.runlocal)
    else:
        run_cmd(' '.join(cmd))

@graphviz(label_prefix="Filtering variants\n", **style_normal)
@transform(ngsEasy_regionFilter, suffix('.targets.vcf'), '.variants.vcf')
@timejob(logger)
def ngsEasy_variantSites(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    cmd = [
        'awk', '-f', os.path.dirname(os.path.realpath(__file__))+'/scripts/rmRefCalls.awk',
        '<', input_file,
        '>', output_file
        ]
    # run job
    if sge:
        run_sge(' '.join(cmd),
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=1, mem=2,runlocally=options.runlocal)
    else:
        run_cmd(' '.join(cmd))

## make a soft filter (PASS for variant that pass)
### make annovar only use the passed variants

#--------------------
# patch to vcf4.0 in-place (required for SAVANT)
#--------------------
@graphviz(label_prefix="Patch VCF\n", **style_normal)
@transform(ngsEasy_variantSites, suffix('.variants.vcf'), '.patched.vcf')
@timejob(logger)
def ngsEasy_patchVcf(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    # convert to bcf
    cmd = [
        'sed',
        '\'s/Number=[A|G|R]/Number=./\'',
        input_file,
        '>', output_file ]
    # run job
    if sge:
        run_sge(' '.join(cmd),
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=1, mem=2,runlocally=options.runlocal)
    else:
        run_cmd(' '.join(cmd))

#--------------------
# Merge filtered variants with GATK or own script
#--------------------
@graphviz(label_prefix="Compress and Index variants\n", **style_normal)
@transform(ngsEasy_patchVcf, regex(r'(.+)\.vcf'), r'\1.vcf.gz.tbi')
@timejob(logger)
def ngsEasy_compressIndexVCF_2(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    # quality filtering
    cmds = [
        [ pipeconfig['software']['bgzip'], '-c', '-f', input_file, '>', input_file+'.gz' ],
        [ pipeconfig['software']['tabix'], input_file+'.gz' ]
            ]
    # run job
    for i, cmd in enumerate(cmds):
        if sge:
            run_sge(' '.join(cmd),
                jobname="_".join([inspect.stack()[0][3], str(i+1), p.RG('SM')]),
                fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
                cpu=1, mem=2,runlocally=options.runlocal)
        else:
            run_cmd(' '.join(cmd))

@follows(ngsEasy_compressIndexVCF_2)
@graphviz(label_prefix="Merging VCF\n", **style_normal)
@collate(ngsEasy_patchVcf, regex(r'(.+)\.(..)\.patched\.vcf'), r'\1.merged.vcf')
@timejob(logger)
def ngsEasy_mergeVcf(input_files,output_file):
    p = loadConfiguration(input_files[0][:input_files[0].rfind('/')])
    # at least supported by 2 callers
    minEvidence = min([2,len(input_files)])
    if pipeconfig['switch']['gatk']:
        cmd = [
            pipeconfig['software']['java'],
            '-Xmx'+str(pipeconfig['resources']['gatk']['CV']['java_mem'])+'g',
            '-Djava.io.tmpdir='+'/'.join([pipeconfig['path']['analysis'], p.wd(), 'tmp']),
            '-jar', pipeconfig['software']['gatk'],
            '-T', 'CombineVariants',
            '-R', pipeconfig['reference']['genome']['sequence'],
            '-genotypeMergeOptions', 'UNIQUIFY',
            #'-genotypeMergeOptions', 'PRIORITIZE',
            #'-genotypeMergeOptions', 'UNSORTED',
            '-minimalVCF',
            '--minimumN', str(minEvidence),
            '-o', output_file ] + \
            [ '--variant:'+input_file[-14:-12]+' '+input_file for input_file in input_files ]
    else:
        cmd = [
            os.path.dirname(os.path.realpath(__file__))+'/scripts/mergevcf.py',
            '-s', pipeconfig['reference']['genome']['dict'],
            '-o', output_file,
            '-e', str(minEvidence) ] + [infile+'.gz' for infile in input_files]
    # run job
    if sge:
        run_sge(' '.join(cmd),
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=pipeconfig['resources']['gatk']['CV']['cpu'], mem=pipeconfig['resources']['gatk']['CV']['mem'],runlocally=options.runlocal)
    else:
        run_cmd(' '.join(cmd))
    # cleanup
    if options.cleanup:
        for infile in input_files:
            zeroFile(infile)

#--------------------
# annotate w/SAVANT and ANNOVAR
#--------------------
@active_if(pipeconfig['switch']['savant'])
@graphviz(label_prefix="Annotate Variants (SAVANT)\n", **style_normal)
@transform([ngsEasy_identifyVariants,ngsEasy_mergeVcf], suffix('.vcf'), '.savant.vcf')
@timejob(logger)
def ngsEasy_savant(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    # MAKE SURE: ALT[1]!="."
    cmd = ' '.join([
        pipeconfig['software']['savant'],
        '-c', pipeconfig['reference']['savant']['config'],
        '-i', input_file,
        '-o', output_file[:output_file.rfind('.')]
        ])
    # run job
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=1, mem=2,runlocally=options.runlocal)
    else:
        run_cmd(cmd)

@active_if(pipeconfig['switch']['snpEff'])
@graphviz(label_prefix="Annotate Variants (snpEff)\n", **style_normal)
@transform([ngsEasy_identifyVariants,ngsEasy_mergeVcf], suffix('.vcf'), '.snpEff.html')
@timejob(logger)
def ngsEasy_snpEff(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    cmd = ' '.join([
        pipeconfig['software']['java'],
        '-Xmx'+str(pipeconfig['resources']['snpEff']['mem'])+'g',
        '-Djava.io.tmpdir='+'/'.join([pipeconfig['path']['analysis'], p.wd(), 'tmp']),
        '-jar', pipeconfig['software']['snpEff'],
        '-dataDir', pipeconfig['reference']['snpEff']['dbdir'],
        '-stats',
        pipeconfig['reference']['snpEff']['database'],
        input_file,
        '>',
        output_file,
        ])
    # run job
    if sge:
        run_sge(cmd,
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=pipeconfig['resources']['snpEff']['cpu'], mem=pipeconfig['resources']['snpEff']['mem'],runlocally=options.runlocal)
    else:
        run_cmd(cmd)


@active_if(pipeconfig['switch']['annovar'])
@graphviz(label_prefix="Annotate Variants (ANNOVAR)\n", **style_normal)
@transform([ngsEasy_identifyVariants,ngsEasy_mergeVcf], suffix('.vcf'), '.annovar.'+ pipeconfig['reference']['annovar']['buildver'] + '_multianno.txt')
@timejob(logger)
def ngsEasy_annovar(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    cmds = []
    # prepare input file
    avinput_file = input_file+'.avinput'
    cmds.append(' '.join([
        pipeconfig['software']['annovar']+'convert2annovar.pl',
        '--includeinfo',
        '--withzyg',
        '--format', 'vcf4',
        input_file,
        '--outfile', avinput_file
        ]))
    # run table_annovar
    _protocol = ','.join(
        pipeconfig['reference']['annovar']['databases']['geneanno'] +
        pipeconfig['reference']['annovar']['databases']['regionanno'] +
        pipeconfig['reference']['annovar']['databases']['filter']
        )
    _operation = ','.join(
        ['g' for x in pipeconfig['reference']['annovar']['databases']['geneanno'] ] +
        ['r' for x in pipeconfig['reference']['annovar']['databases']['regionanno'] ] +
        ['f' for x in pipeconfig['reference']['annovar']['databases']['filter'] ]
        )
    cmds.append(' '.join([
        pipeconfig['software']['annovar']+'table_annovar.pl',
        '--protocol', _protocol,
        '--operation', _operation,
        #'--csvout',
        #'--otherinfo',
        '--buildver', pipeconfig['reference']['annovar']['buildver'],
        '--remove',
        '--outfile', output_file[:output_file.find('.'+pipeconfig['reference']['annovar']['buildver']+'_multianno')],
        avinput_file,
        pipeconfig['reference']['annovar']['directory']
        ]))
    # run job
    for cmd in cmds:
        if sge:
            run_sge(cmd,
                jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
                fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
                cpu=1, mem=2,runlocally=options.runlocal)
        else:
            run_cmd(cmd)

#--------------------
# coverage
#--------------------
@graphviz(label_prefix="Base coverage\n", **style_normal)
@transform(ngsEasy_markDuplicates, suffix('.dupemk.bam'), '.coverageBed')
@timejob(logger)
def ngsEasy_coverageBed(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    # reduce to target regions of targeted sequencing or exome
    cmd = [
        pipeconfig['software']['bedtools'],
        'coverage',
        '-d',
        '-abam', input_file,
        '-b', ngsconfig['ngsanalysis'][p.ngsAnalysis]['tas']['covered_bed'],
        '>', output_file
        ]
    # run job
    if sge:
        run_sge(' '.join(cmd),
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=1, mem=2,runlocally=options.runlocal)
    else:
        run_cmd(' '.join(cmd))

@active_if(pipeconfig['switch']['gatk'])
@graphviz(label_prefix="Covered intervals\n", **style_normal)
@transform(ngsEasy_markDuplicates, suffix('.dupemk.bam'), '.coveredIntervals')
@timejob(logger)
def ngsEasy_coveredIntervals(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    # reduce to target regions of targeted sequencing or exome
    cmd = [
        pipeconfig['software']['java'],
        '-Xmx'+str(pipeconfig['resources']['gatk']['FCI']['java_mem'])+'g',
        '-Djava.io.tmpdir='+'/'.join([pipeconfig['path']['analysis'], p.wd(), 'tmp']),
        '-jar', pipeconfig['software']['gatk'],
        '-T', 'FindCoveredIntervals',
        '-R', pipeconfig['reference']['genome']['sequence'],
        '-I', input_file,
        '-o', output_file,
        '--coverage_threshold', str(4)
        ]
    if sge:
        run_sge(' '.join(cmd),
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=pipeconfig['resources']['gatk']['FCI']['cpu'], mem=pipeconfig['resources']['gatk']['FCI']['mem'],runlocally=options.runlocal)
    else:
        run_cmd(' '.join(cmd))

'''
__FUTURE__
#run PCA-SVM for coverage variant detection (per gene?)
'''

#--------------------
# Structural variation detection (mapping based)
#--------------------

@active_if(pipeconfig['switch']['SV'])
@graphviz(label_prefix="SV delly\n", **style_normal)
@transform(ngsEasy_markDuplicates, suffix('.dupemk.bam'), '.SV.delly')
@timejob(logger)
def ngsEasy_delly(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    cmd = [
        # pipeconfig['software']['java'],
        # '-Xmx'+str(pipeconfig['resources']['gatk']['FCI']['java_mem'])+'g',
        # '-Djava.io.tmpdir='+'/'.join([pipeconfig['path']['analysis'], p.wd(), 'tmp']),
        # '-jar', pipeconfig['software']['gatk'],
        # '-T', 'FindCoveredIntervals',
        # '-R', pipeconfig['reference']['genome']['sequence'],
        # '-I', input_file,
        # '-o', output_file,
        # '--coverage_threshold', str(4)
        ]
    if sge:
        run_sge(' '.join(cmd),
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=pipeconfig['resources']['delly']['cpu'], mem=pipeconfig['resources']['delly']['mem'],runlocally=options.runlocal)
    else:
        run_cmd(' '.join(cmd))


@active_if(pipeconfig['switch']['SV'])
@graphviz(label_prefix="SV exomedepth\n", **style_normal)
@transform(ngsEasy_markDuplicates, suffix('.dupemk.bam'), '.SV.exomedepth')
@timejob(logger)
def ngsEasy_exomedepth(input_file,output_file):
    p = loadConfiguration(os.path.abspath(os.path.dirname(input_file)))
    cmd = [
        pipeconfig['software']['R'],
        os.path.dirname(os.path.realpath(__file__))+'/scripts/exomedepth.R',
        input_file,
        output_file,
        ]
    if sge:
        run_sge(' '.join(cmd),
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=pipeconfig['resources']['pindel']['cpu'], mem=pipeconfig['resources']['pindel']['mem'],runlocally=options.runlocal)
    else:
        run_cmd(' '.join(cmd))


@active_if(pipeconfig['switch']['SV'])
@graphviz(label_prefix="SV exon homozygosity\n", **style_normal)
@transform(ngsEasy_mergeVcf, suffix('.dupemk.bam'), '.SV.exonhomo.bed')
@timejob(logger)
def ngsEasy_exonhomo(input_file,output_file):
    p = loadConfiguration(os.path.abspath(os.path.dirname(input_file)))
    cmd = [
        pipeconfig['software']['R'],
        os.path.dirname(os.path.realpath(__file__))+'/scripts/exonHomo.py',
        input_file,
        ngsconfig['ngsanalysis'][p.ngsAnalysis]['tas']['covered_bed'],
        output_file
        ]
    if sge:
        run_sge(' '.join(cmd),
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=1, mem=2,runlocally=options.runlocal)
    else:
        run_cmd(' '.join(cmd))

@active_if(pipeconfig['switch']['SV'])
@graphviz(label_prefix="SV exon homozygosity\n", **style_normal)
@transform(ngsEasy_alignment, suffix('.dupemk.bam'), '.SV.exonhomo.bed')
@timejob(logger)
def ngsEasy_SLOPE(input_file,output_file):
    p = loadConfiguration(os.path.abspath(os.path.dirname(input_file)))
    cmd = [
        # http://www-genepi.med.utah.edu/suppl/SLOPE/index.html
        ]
    if sge:
        run_sge(' '.join(cmd),
            jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
            fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
            cpu=1, mem=2,runlocally=options.runlocal)
    else:
        run_cmd(' '.join(cmd))
#
### breakdancer?
##### CNVnator
####### m-HMM
##### http://www-genepi.med.utah.edu/suppl/SLOPE/index.html
###
#

'''Summarize large scale structural variants'''
@timejob(logger)
@graphviz(label_prefix="Summarize SV\n", **style_collect)
@collate([ngsEasy_delly, ngsEasy_exomedepth, ngsEasy_exonhomo, ngsEasy_SLOPE],
    regex(r'(.+\.SV)\.\S+'), r'\1.summary')
@timejob(logger)
def ngsEasy_summarizeSV(input_files,output_file):
    with open(output_file,'w') as outfh:
        for infile in flatten(input_files):
            if infile.endswith('pdf') or infile.endswith('coverage'):
                continue
            with open(infile) as infh:
                fields, ordering, postfix = parsePicard(infh)
                # write summary to file
                for n in sorted(ordering.keys()):
                    k = ordering[n]
                    v = fields[k]
                    if len(v) == 0:
                        continue
                    for i,p in enumerate(postfix):
                        try:
                            float(v[i])
                        except ValueError:
                            pass
                        else:
                            print >> outfh, '{:<30} {:<16} {:<26} {}'.format(infile[infile.rfind('.')+1:], postfix[i], k, v[i])

#--------------------
# final target
#--------------------
@graphviz(label_prefix="END\n", **style_end)
@merge([ngsEasy_readQC,
    ngsEasy_analyzeCovariates,
    ngsEasy_alignmentQC,
    ngsEasy_coverageBed,
    ngsEasy_coveredIntervals,
    ngsEasy_varStats,
    ngsEasy_savant,
    ngsEasy_annovar,
    ngsEasy_snpEff,
    ngsEasy_summarizeSV],
    'final.checkpoint')
@timejob(logger)
def ngsEasy_final(input_files, output_file):
    # just touch the out_file
    with open(output_file, 'w') as fh:
        for infile in input_files:
            if os.stat(infile).st_size > 0:
                print >> fh, infile
    logger.info('##############################')
    logger.info('########## HEUREKA! ##########')
    logger.info('##############################')


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   setup target/forced tasks
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
if options.forcedtasks:
    try:
        forcedtasks = [ globals()[x.strip()] for x in options.forcedtasks.split(',')]
    except:
        print >> sys.stderr, "valid tasks are:\n\t",'\n\t'.join([ x for x in locals().keys() if x.startswith('ngsEasy_') ])
        logger.critical("unkown forced task in %s" % options.forcedtasks)
        sys.exit(1)
    logger.info("FORCED_TASKS: "+','.join([ getattr(f,'func_name') for f in forcedtasks]))
else:
    forcedtasks = []

if options.targettasks:
    try:
        targettasks = [ globals()[x.strip()] for x in options.targettasks.split(',')]
    except:
        print >> sys.stderr, "valid tasks are:\n\t",'\n\t'.join([ x for x in locals().keys() if x.startswith('ngsEasy_') ])
        logger.critical("unkown target task in %s" % options.targettasks)
        sys.exit(1)
    logger.info("TARGET_TASKS: "+','.join([ getattr(f,'func_name') for f in targettasks]))
else:
    targettasks = []

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Print list of tasks
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
if options.just_print:
    pipeline_printout(sys.stdout, [], verbose=options.verbose, verbose_abbreviated_path = 3)

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Print flowchart
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
elif options.flowchart:
    # use file extension for output format
    output_format = os.path.splitext(options.flowchart)[1][1:]
    pipeline_printout_graph(open(options.flowchart, "w"), output_format, [ngsEasy_final],
        no_key_legend = False, minimal_key_legend=True, draw_vertically=True,
        pipeline_name='ngsEasy - ngs pipeline', user_colour_scheme={"colour_scheme_index" :3,
                                                "Pipeline"      :{"fontcolor" : '"#FF3232"' },
                                                "Key"           :{"fontcolor" : "Red",
                                                                  "fillcolor" : '"#F6F4F4"' },
                                                "Task to run"   :{"linecolor" : '"#0044A0"' },
                                                "Final target"  :{"fillcolor" : '"#EFA03B"',
                                                                  "fontcolor" : "black",
                                                                  "dashed"    : 0           },
                                                "Up-to-date task"  :{"fillcolor" : "yellow",
                                                                  "fontcolor" : "black",
                                                                  "dashed"    : 1           }})

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Run Pipeline
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
else:
    pipeline_run(target_tasks = targettasks, forcedtorun_tasks = forcedtasks,
        multithread = options.jobs if sge else 1, multiprocess = 1 if sge else options.jobs,
        logger = logger, verbose=options.verbose, verbose_abbreviated_path = 3,
        checksum_level = 1 if options.debug else 3, touch_files_only=options.touch)

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Cleanly end drmaa session
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
if sge:
    sge.exit()


'''
VARIANT RECALIBRATOR - ONLY WORKS IF ENOUGH VARIANTS PRESENT
@transform(ngsEasy_haplotypeCaller,suffix('.vcf'),'.recalibrated.vcf')
def ngsEasy_variantRecalibratorHC(input_file, output_file):
    basepath = input_file[:input_file.rfind('/')]
    p = loadConfiguration(basepatg)
    recal_file = basepath+'/'+input_file+'.recal'
    tranches_file = basepath+'/'+input_file+'.tranches'
    # recalibrateSNP, apply, recalibrateIndels, apply
    cmd = [
        " ".join([
            pipeconfig['software']['java'],
            '-Xmx'+str(pipeconfig['resources']['gatk']['VR']['java_mem'])+'g',
            '-Djava.io.tmpdir='+'/'.join([pipeconfig['path']['analysis'], p.wd(), 'tmp']),
            '-jar', pipeconfig['software']['gatk'],
            '-T', 'VariantRecalibrator',
            '-R', pipeconfig['reference']['genome']['sequence'],
            '-nt', str(pipeconfig['resources']['gatk']['VR']['nt']),
            '-input', input_file,
            '-resource:dbsnp,known=true,training=false,truth=false,prior=6.0', pipeconfig['reference']['snps'][0],
            '-resource:hapmap,known=false,training=true,truth=true,prior=15.0', pipeconfig['reference']['snps'][1],
            '-resource:omni,known=false,training=true,truth=false,prior=12.0', pipeconfig['reference']['snps'][2],
            '-an QD', '-an HaplotypeScore', '-an MQRankSum', '-an ReadPosRankSum', '-an FS', '-an MQ', '-an InbreedingCoeff',
            '-mode', 'BOTH',
            '-recalFile', recal_file,
            '-tranchesFile', tranches_file,
            '-rscriptFile', +input_file+'.recalplots.R'
        ]),
        " ".join([
            pipeconfig['software']['java'],
            '-Xmx'+str(pipeconfig['resources']['gatk']['VR']['java_mem'])+'g',
            '-Djava.io.tmpdir='+'/'.join([pipeconfig['path']['analysis'], p.wd(), 'tmp']),
            '-jar', pipeconfig['software']['gatk'],
            '-T', 'ApplyRecalibration',
            '-R', pipeconfig['reference']['genome']['sequence'],
            '-input', input_file,
            '--ts_filter_level', '99.0',
           '-tranchesFile', tranches_file,
           '-recalFile', recal_file,
           '-mode', 'BOTH',
           '-o', output_file
        ])
    ]
    # run recalibration
    for cmd in cmds:
        if sge:
            run_sge(cmd,
                jobname="_".join([inspect.stack()[0][3], p.RG('SM')]),
                fullwd=pipeconfig['path']['analysis']+'/'+p.wd(),
                cpu=1, mem=2)
        else:
            run_cmd(cmd)
#@transform(ngsEasy_variantRecalibratorHC, suffix('.vcf'), '.filtered.vcf')

VCFTOOLS VARIANT FILTER - NOW USING COLLATED FILES WITH BCFTOOLS FILTER
@transform(ngsEasy_haplotypeCaller, suffix('.vcf'), '.filtered.vcf')
def ngsEasy_siteFilterHC(input_file,output_file):
    p = loadConfiguration(input_file[:input_file.rfind('/')])
    cmd = " ".join([
        pipeconfig['software']['vcftools'],
        '--vcf',input_file,
        '--recode',
        '--recode-INFO-all',
        '--min-meanDP', str(ngsconfig['ngsanalysis'][p.ngsAnalysis]['vcf']['minDP']),
        '--minQ', str(ngsconfig['ngsanalysis'][p.ngsAnalysis]['vcf']['minQ']),
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



