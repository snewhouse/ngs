#!/usr/bin/env python

__doc__=='''
testing ground for ruffus
'''

import sys
from ruffus import *
import time


inputfiles = [ './sand/coco.fastq', './sand/caca.fastq' ]

sampledir = 'sand'

##### PIPELINE FUNCTIONS
##### PIPELINE FUNCTIONS
##### PIPELINE FUNCTIONS
##### PIPELINE FUNCTIONS

@mkdir(sampledir)
@originate(inputfiles)
def startingFiles(output_file):
    print 'ORIGINATE'
    print output_file
    # create both files as necessary
    with open(output_file, "w") as oo: pass
    return

@transform(startingFiles, suffix('.fastq'), '.sam')
def leftTask(input_file,output_file):
    print 'LEFT'
    print input_file
    print output_file
    with open(output_file,'w') as fh:
        fh.write('dummy')
        pass


@transform(startingFiles,suffix('.fastq'),'.sue')
def riteTask(input_file,output_file):
    print 'RITE'
    print input_file
    print output_file
    with open(output_file,'w') as fh:
        fh.write('dummy')
        pass

#@transform([leftTask,riteTask], regex(r'.*?([^/]+)\.(s..)'),[r'\1.\2.Success'])
#@merge([leftTask,riteTask], 'output.file')
@collate([leftTask,riteTask], regex(r'.*?([^/]+)\.(s..)'),[r'\1.Success'])
def testMe(input_files, output_files):
    print type(input_files)
    print '### IN  ###', input_files
    print '### OUT ###', output_files

# run pipeline
pipeline_run()
#pipeline_printout()



