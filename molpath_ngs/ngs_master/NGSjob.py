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

 Concat ReadGroup_id_RGID
    ReadGroup_library_RGLB
 pull from fastq ReadGroup_platform_RGPL
 pull from fastq ReadGroup_platform_unit_RGPU
     ReadGroup_SeqCentre_RGCN
     ReadGroup_Desc_RGDS
     ReadGroup_runDate_RGDT
XXX QUAL
ANALYSIS Pipeline
XXX PE
bed_list
bed_type
XXX email_bioinf
PIPECONFIG Fastq_dir
PIPECONFIG BAM_dir
PIPECONFIG pipeline_dir
XXX sleep
'''

class Patients(list):
    def __init__(self, fh, delim='\t', strict=True):
        self.fields = []  # must be concordant to patient slots
        # parser with sanity check
        for lineNum, line in enumerate(fh):
            if lineNum == 0:
                # header
                self.fields = line.rstrip().split(delim)
                if strict:
                    #check if header equal slots
                    try:
                        assert all([ x[0]==x[1] for x in zip(self.fields, Patient.__slots__) ])
                    except AssertionError:
                        logging.error("Fields in configuration file are wrong")
            else:
                self.append(Patient(line))
        return


class Patient:

    __slots__ = ("Fastq_file_prefix","ReadGroup_sample_RGSM",
        "ReadGroup_id_RGID","ReadGroup_library_RGLB",
        "ReadGroup_platform_RGPL","ReadGroup_platform_unit_RGPU",
        "ReadGroup_SeqCentre_RGCN","ReadGroup_Desc_RGDS",
        "ReadGroup_runDate_RGDT","QUAL","Pipeline","PE","bed_list",
        "bed_type","email_bioinf","Fastq_dir","BAM_dir","pipeline_dir",
        "sleep"
        # additional attributes
        )

    def __init__(self, sline):
        if isinstance(sline, dict):
            try:
                assert len(sline) == len(self.__slots__)
            except:
                print >> sys.stderr, '##', sline, '##'
                raise Exception('Patient data parsing error (from JSON)')
            # read fields
            for k,v in sline.items():
                setattr(self, k, v)
        else:
            args = sline.strip().split("\t")
            try:
                assert len(args) == len(self.__slots__)
            except:
                print >> sys.stderr, '##', sline, '##'
                raise Exception('Patient data parsing error')
            # read fields
            for i in range(len(args)):
                setattr(self, Patient.__slots__[i], args[i])
        return

    def __repr__(self):
        fields = zip(Patient.__slots__, \
         [ getattr(self, Patient.__slots__[i]) for i in range(len(Patient.__slots__)) ])
        #print >> sys.stderr, self.wd()
        return "\n".join(['{1:>30} {0:>2} {2:<40}'.format(i, *f) for i,f in enumerate(fields)])

    def __str__(self):
        return '< Patient: {:_^50} {:<30} {:<30} {:<10} >'.format(
            getattr(self, Patient.__slots__[2]),
            getattr(self, Patient.__slots__[3]),
            getattr(self, Patient.__slots__[1]),
            getattr(self, Patient.__slots__[8]),
            )

    def __getitem__(self, key):
        return getattr(self, key)

    def __lt__(self, other):  # sort by time year-month-day
        try:
            selfDate = tuple(map(int, self.ReadGroup_runDate_RGDT.split("-")[:3]))
            otherDate = tuple(map(int, other.ReadGroup_runDate_RGDT.split("-")[:3]))
        except:
            logging.warning("date sorting error: not YYYY-MM-DD")
            return True
        return selfDate < otherDate

    def save(self,outfile):
        with open(outfile, 'w') as outf:
            json.dump(self.__dict__, outf)
        return

    # workdir (BAM_)
    def wd(self):
        return '/'.join([self.BAM_dir.rstrip('/'), self.uniqueID()])

    def fastq(self):
        return [ '/'.join([self.Fastq_dir,x.strip()]) for x in self.Fastq_file_prefix.split(',') ]

    def linkedfastq(self):
        linkedfq = []
        for i, x in enumerate(self.fastq()):
            linkedfq.append('/'.join([self.wd(), self.RGSM()+'_'+str(i+1)+'.fastq']))
        return linkedfq

    def RGSM(self):
        return self.ReadGroup_sample_RGSM

    def RGID(self):
        return self.ReadGroup_id_RGID

    def RGLB(self):
        return self.ReadGroup_library_RGLB

    def RGPL(self):
        return self.ReadGroup_platform_RGPL.lower()

    def RGPU(self, setnew=None):
        if setnew:
            self.ReadGroup_platform_unit_RGPU = setnew
        return self.ReadGroup_platform_unit_RGPU

    def RGCN(self):
        return self.ReadGroup_SeqCentre_RGCN

    def RGDS(self):
        return self.ReadGroup_Desc_RGDS

    def RGDT(self):
        return self.ReadGroup_runDate_RGDT

    def uniqueID(self):  # sample_name
        return self.RGSM() + self.RGLB() + self.RGPL() + self.RGDT()

    def isPE(self):
        try:
            self.PE = int(self.PE)
            assert self.PE in [0,1]
        except AssertionError:
            raise Exception("PE parameter not 0 or 1")
        except:
            raise

        if self.PE == 1:
            return True
        return False

    def librarytype(self):
        try:
            self.PE = int(self.PE)
            assert self.PE in [0,1]
        except AssertionError:
            raise Exception("PE parameter not 0 or 1")
        except:
            raise

        if self.PE == 1:
            return 'PE'
        return 'SE'

    def analysis(self):
        return 'default'


if __name__ == "__main__":
    patients = Patients(sys.stdin)
    for patient in patients:
        patient.save('sample.json')
        with open('sample.json') as fh:
            js = json.load(fh)
            patient2 = Patient(js)
        print patient
        print patient2
        print repr(patient)
