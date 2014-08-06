#!/usr/bin/env python

__doc__=="""
#############################################################################################
# NGS pipeline for molecular diagnostics (can be dockerized)                                #
# -- Organisation: Viapath LLC (and KCL/SLaM/NHS)                                           #
# -- Date: 06/08/2014                                                                       #
#############################################################################################
"""
__author__ = "David Brawand, Stephen Newhouse, Amos Folarin, Aditi Gulati"
__copyright__ = "LGPL"
__credits__ = ['Stephen Newhouse', 'Amos Folarin', 'Aditi Gulati']
__license__ = "LGPL"
__version__ = "1.4"
__maintainer__ = "David Brawand"
__email__ = "dbrawand@nhs.net, stephen.j.newhouse@gmail.com, amosfolarin@gmail.com, aditi.gulati@nhs.net"
__status__ = "Development"  # ["Prototype", "Development",  "Production"]

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   basdic imports and paths
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

import os, sys, re
import subprocess

exe_path = os.path.split(os.path.abspath(sys.argv[0]))[0]
#sys.path.insert(0,os.path.abspath(os.path.join(exe_path, "..", "..")))

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   options and input parsing
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
from optparse import OptionParser

usage = "usage: %prog [options] <sample descriptions>"
parser = OptionParser(version="%prog 0.1", usage=usage)

#   mandatory options
parser.add_option("-p", dest="pipepar", default='ngs_config.json',metavar="STRING", help="Pipeline parameters (ngs_config.json)")
parser.add_option("-l", dest="logfile", default=None,             metavar="STRING", help="logfile (default STDERR)")
#   general options: verbosity / logging
parser.add_option("-v", dest="verbose", action="count", default=0, help="Print more detailed messages (eg. -vvv)")
#   pipeline
parser.add_option("-j", "--jobs", dest="jobs",
                  default=1,
                  metavar="jobs",
                  type="int",
                  help="Specifies the number of jobs (operations) to run in parallel.")
parser.add_option("--flowchart", dest="flowchart",
                  metavar="FILE",
                  type="string",
                  help="Print flowchart of the pipeline to FILE. Flowchart format "
                       "depends on extension. Alternatives include ('.dot', '.jpg', "
                       "'*.svg', '*.png' etc). Formats other than '.dot' require "
                       "the dot program to be installed (http://www.graphviz.org/).")
parser.add_option("-n", "--just_print", dest="just_print",
                    action="store_true", default=False,
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
from Patient import Patients

# read pipeline parameters
logger.info("reading configuration from %s", options.pipepar)
with open(options.pipepar) as pp:
    pipeconfig = json.load(pp)

# read patient files
for patientfile in args:
    # read sample/patient data
    with open(patientfile) as fh:
        patients = Patients(fh)

    # print patient data in readable format
    for i, patient in enumerate(patients):
        print i, patient
        print repr(patient)

        # make target directories and write configuration
        try:
            os.mkdir(patient.wd())
        except OSError:
            pass
        else:
            # write configurations
            configfile = patient.wd()+'/.molpath'
            patient.save(configfile)

        # check fastQ files and symlink to sample file

        # set input files


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
        patient_config = Patient(json.load(fh))
    return patient_config

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

# this filters the fastq
@transform(inputfiles,
            formatter("(?P<PAIR>[^\.]+)\.fastq", "(?P<PAIR>[^\.]+)\.fastq"),
            ["{PAIR[0]}.filtered.fastq", "{PAIR[0]}.unpaired.fastq",
             "{PAIR[1]}.filtered.fastq", "{PAIR[1]}.unpaired.fastq"])
def trimmomatic(input_files, output_files):
    # load configurations
    p = loadConfiguration(input_files[0])
    # configure job
    job_name = 'trimmomatic'+p.uid()
    cmd = ' '.join([
        config['software']['java_v1_7'] + '/java',
        '-XX:ParallelGCThreads=4 -Xmx', str(config['resources']['trimmomatic']['java_mem'])+'g',
        '-jar', config['software']['trimmomatic'],
        p.librarytype(),
        '-phred64',
        '-trimlog', p.uid()+'.trimming.log',
        '-threads', config['resources']['trimmomatic']['cpu'],
        ' '.join(input_files),
        ' '.join(output_files),
        ' '.join(config['pipeline'][p.analysis()]['trimmomatic']) ])
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
        print cmd

    for ofi in output_files:
        with open(ofi,'w') as fh:
            fh.write('dummy')
            pass
    print input_files, output_files
    return

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
    pipeline_printout(sys.stdout, [combineBlastResults], verbose=options.verbose)


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Print flowchart
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
elif options.flowchart:
    # use file extension for output format
    output_format = os.path.splitext(options.flowchart)[1][1:]
    pipeline_printout_graph (open(options.flowchart, "w"),
                             output_format,
                             [combineBlastResults],
                             no_key_legend = True)
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Run Pipeline
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
else:
    pipeline_run([combineBlastResults],  multiprocess = options.jobs,
                        logger = logger, verbose=options.verbose)

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Cleanly end drmaa session
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
if sge:
    sge.exit()


