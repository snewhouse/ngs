#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import sys
import datetime
import json
import hashlib

'''
ProjectID
SampleID/PatientID (RGSM)
FASTQ1
FASTQ2

NGStype     (RGPL, WGS/WEX/TGS, RGLB)
NGSanalysis (BED, analysis)

ReadGroup_SeqCentre_RGCN
ReadGroup_Desc_RGDS
ReadGroup_runDate_RGDT

AUTOMATIC RGPU (from fastq)
AUTOMATIC RGID (RGSM.NGStype.NGSanalysis.${DATE})
AUTOMATIC PE

IMPLICIT RGLB
IMPLICIT RGPL
IMPLICIT RGDS
'''

# parses NGSjobs file (checks field order)
class NGSjobs(list):

    __aliases__ = {
        'runDate': 'RGDT',
        'sampleID': 'RGSM',
        'worksheetID': 'runID',
        'samplesheet': 'sampleSheet'
     } # field aliases

    def __init__(self, fh, delim='\t'):
        self.fields = []  # must be concordant to patient slots
        # parser with sanity check
        for lineNum, line in enumerate(fh):
            if lineNum == 0:
                # header replace aliases
                self.fields = [ NGSjobs.__aliases__[x] if x in NGSjobs.__aliases__.keys() \
                    else x for x in line.rstrip().split(delim) ]
            else:
                # deserialize into hash
                data = dict(zip(self.fields, line.rstrip().split(delim)))
                self.append(NGSjob(data))
        return


class NGSjob:

    __slots__ = (
        "projectID",  #dev/diag/res
        "RGSM",  #sampleID
        "runID",  # worksheetID
        "sampleSheet", # runconfig
        "ngsType",  # RGPL(ILLUMINA), RGLB(WGS/WEX/TGS), RGDS(bait_liver), RGCN(molpath)
        "ngsAnalysis",  # (BED, analysis)
        "cleanup",  # remove intermediary files
        "FASTQ1",  # from sampleSheet
        "FASTQ2",  # from sampleSheet
        "RGPU",  # from FASTQ
        "RGCN",  # from ngsType
        "RGLB",  # from ngsType
        "RGPL",  # from ngsType
        "RGDS",  # from ngsType
        "RGDT",  # from sampleSheet
        "RGID",  # from (RGSM.NGStype.NGSanalysis.RGDT)
        "_inputfiles",  # the FASTQ files to start with, can be anonymized filename
        "_sample",  # from filename
        "_run_id",  # from filename
        "_lane",  # from filename
        "_pair",  # from filename
        "githash",  # git revision history
        "_pipeconfig"  # pipeline configration (reference and software at least)
        )

    __required__ = (
        "projectID",  #dev/diag/res
        "RGSM",  #sampleID
        "runID",  # worksheetID
        "RGDT", # from sampleSheet
        "RGPU", # from FASTQ
        "RGCN", # from ngsType
        "RGLB", # from ngsType
        "RGPL", # from ngsType
        "RGDS", # from ngsType
        #"ngsType",  # RGPL(ILLUMINA), RGLB(WGS/WEX/TGS), RGDS(bait_liver), RGCN(molpath)
        "ngsAnalysis",  # (BED, analysis)
        "FASTQ1"  # from sampleSheet
    )

    # GENERIC INIT
    def __init__(self, *data):
        for i, d in enumerate(data):
            # check if keys are valid
            try:
                unkown = set(d.keys()).difference(set(self.__slots__))
                assert len(unkown)==0
            except AssertionError:
                    raise Exception('ERROR: configuration field unkown: ' + ",".join(unkown))
            # set initial values (None)
            for s in self.__slots__:
                try:
                    setattr(self, s, d[s])
                except KeyError:
                    if i==0:
                        setattr(self, s, None)
                    else:
                        pass
        return

    def __repr__(self):
        fields = zip(NGSjob.__slots__, \
         [ getattr(self, NGSjob.__slots__[i]) for i in range(len(NGSjob.__slots__)) ])
        #print >> sys.stderr, self.wd()
        return "\n".join(['{1:>30} {0:>2} {2:<40}'.format(i, *f) for i,f in enumerate(fields)]+[''])

    def __str__(self):
        try:
            return '< NGSjob: {:}:{:}:{:} {:<10} >'.format(
                getattr(self, NGSjob.__required__[0]),
                getattr(self, NGSjob.__required__[1]),
                getattr(self, NGSjob.__required__[2]),
                getattr(self, NGSjob.__required__[3]),
                )
        except:
            print >> sys.stderr, dir(self)
            raise

    def __getitem__(self, key):
        return getattr(self, key)

    def save(self,outfile):  # write to JSON
        with open(outfile, 'w') as outf:
            json.dump(self.__dict__, outf)
        return

    # Dependent on slot names
    def __lt__(self, other):  # sort by time year-month-day
        try:
            selfDate = tuple(map(int, self.runDate.split("-")[:3]))
            otherDate = tuple(map(int, other.runDate.split("-")[:3]))
        except:
            return True
        return selfDate < otherDate

    def wd(self):  # work dir
        return self.projectID+'/'+self.RGSM + '/' + self.runID

    def setRGID(self,rgid=None):  # unique identifyer
        if rgid:
            self.RGID = rgid
        elif not self.RGID:
            try:
                assert self.projectID and self.RGSM and self.runID and self.RGDT
            except AssertionError:
                raise Exception("Cannot generate unique identifyer (RGID)")
            except:
                raise
            else:
                self.RGID = hashlib.md5(''.join([self.projectID, self.RGSM, self.runID, self.RGDT])).hexdigest()[:8]
        else:
            pass  # won't overwrite
        return

    def fastq(self,basepath=None):
        if basepath:
            return [ '/'.join([basepath.rstrip('/'), self.FASTQ1]), \
                     '/'.join([basepath.rstrip('/'), self.FASTQ2]) ]
        else:
            return [ self.FASTQ1, self.FASTQ2 ]

    def inputfiles(self,basepath=None):
        if basepath:  # set the inputfiles
            self._inputfiles =  [ '/'.join([basepath.rstrip('/'), self.RG('ID')+"_1.fastq"]), \
                     '/'.join([basepath.rstrip('/'), self.RG('ID')+"_2.fastq"]) ]
        return self._inputfiles

    # get/set ReadGroup
    def RG(self,field):
        if field.upper() == 'SM':
            return self.RGSM
        if field.upper() == 'ID':
            if self.RGID:
                return self.RGID
            return "_".join([self.RGSM, self.ngsType, self.ngsAnalysis, self.RGDT])
        if field.upper() == 'DT':
            return self.RGDT
        if field.upper() == 'PL':
            return self.RGPL
        if field.upper() == 'LB':
            return self.RGLB
        if field.upper() == 'DS':
            return self.RGDS
        if field.upper() == 'PU':
            return self.RGPU
        if field.upper() == 'CN':
            return self.RGCN
        else:
            raise Exception('unknown RG field: %s' % field)

    def librarytype(self):
        if self.FASTQ2:
            return 'PE'
        return 'SE'

    # returns true if requirements are met
    def unmetRequirements(self):
        unmet = []
        for r in NGSjob.__required__:
            if not getattr(self, r):
                unmet.append(r)
        return unmet

if __name__ == "__main__":
    jobs = NGSjobs(sys.stdin)
    for job in jobs:
        job.save('sample.json')
        with open('sample.json') as fh:
            js = json.load(fh)
            job2 = NGSjob(js)
            job3 = NGSjob(js,{"something":"", "sampleID":"field", "RGDS": "new_center", 'RGSM': "different sample"})
        job3.save('sample.json')
        with open('sample.json') as fh:
            js = json.load(fh)
            job4 = NGSjob(js)
        print repr(job)
        print repr(job2)
        print repr(job3)
        print repr(job4)
