#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import sys
import datetime
import json

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
    def __init__(self, fh, delim='\t'):
        self.fields = []  # must be concordant to patient slots
        # parser with sanity check
        for lineNum, line in enumerate(fh):
            if lineNum == 0:
                # header
                self.fields = line.rstrip().split(delim)
            else:
                # deserialize into hash
                data = dict(zip(self.fields, line.rstrip().split(delim)))
                self.append(NGSjob(data))
        return


class NGSjob:

    __slots__ = (
        "projectID",  #dev/diag/res
        "sampleID",  #RGSM
        "runID",  # worksheetID
        "sampleSheet", # runconfig
        "ngsType",  # RGPL(ILLUMINA), RGLB(WGS/WEX/TGS), RGDS(bait_liver), RGCN(molpath)
        "ngsAnalysis",  # (BED, analysis)
        "cleanup",  # remove intermediary files
        "FASTQ1",  # from sampleSheet
        "FASTQ2",  # from sampleSheet
        "RGPU", # from FASTQ
        "RGCN", # from ngsType
        "RGLB", # from ngsType
        "RGPL", # from ngsType
        "RGDS", # from ngsType
        "RGDT", # from sampleSheet
        "RGID" # from (RGSM.NGStype.NGSanalysis.RGDT)
        )

    # GENERIC INIT
    def __init__(self, *data):
        for i, d in enumerate(data):
            # check if keys are valid
            try:
                unkown = set(d.keys()).difference(set(self.__slots__))
                assert len(unkown)==0
            except AssertionError:
                    raise Exception('ERROR: configuration field: ' + ",".join(unkown))
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
        return '< NGSjob: {:_^50} {:<30} {:<30} {:<10} >'.format(
            getattr(self, NGSjob.__slots__[1]),
            getattr(self, NGSjob.__slots__[0]),
            getattr(self, NGSjob.__slots__[4]),
            getattr(self, NGSjob.__slots__[6]),
            )

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
            logging.warning("date sorting error: not YYYY-MM-DD")
            return True
        return selfDate < otherDate

    def wd(self):  # work dir
        return self.projectID+'/'+self.sampleID

    def fastq(self):
        return [ self.FASTQ1, self.FASTQ2 ]

    # get/set ReadGroup
    def RG(self,field):
        if field.upper() == 'SM':
            return self.sampleID
        if field.upper() == 'ID':
            return "_".join([self.sampleID, self.ngsType, self.ngsAnalysis. self.RGDT])
        if field.upper() == 'DT':
            return self.runDate
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


if __name__ == "__main__":
    jobs = NGSjobs(sys.stdin)
    for job in jobs:
        job.save('sample.json')
        with open('sample.json') as fh:
            js = json.load(fh)
            job2 = NGSjob(js)
            job3 = NGSjob(js,{"something":"", "sampleID":"field", "RGDS": "new_center", 'sampleID': "different sample"})
        job3.save('sample.json')
        with open('sample.json') as fh:
            js = json.load(fh)
            job4 = NGSjob(js)
        print repr(job)
        print repr(job2)
        print repr(job3)
        print repr(job4)
