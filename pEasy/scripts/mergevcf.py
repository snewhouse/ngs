#!/usr/bin/env python

__doc__='''
##############################################################
## mergevcf.py | merges tabix indexed vcfs                  ##
## -- a naive and free alternative to GATKs CombineVariants ##
##############################################################
'''
__author__ = "David Brawand"
__credits__ = ['David Brawand']
__license__ = "MIT licence"
__version__ = "0.1"
__maintainer__ = "David Brawand"
__email__ = "dbrawand@nhs.net"
__status__ = "Development"  # ["Prototype", "Development",  "Production"]


import sys
import re
import tabix  # pytabix
import gzip
from collections import Counter

class Format(object):
    def __init__(self,line):
        m = re.match(r'##FORMAT=<ID=([^,]+),Number=([^,]+),Type=([^,]+),Description="([^"]+)">',line)
        try:
            assert m
        except AssertionError:
            raise
        else:
            self.ID = m.group(1)
            self.Number = m.group(2)
            self.Type = m.group(3)
            self.Description = m.group(4)
        return

    def update(self,other):
        if self.ID == other.ID:
            # update Number
            if self.Number == other.Number or (not self.ID.startswith('.')) != (not other.ID.startswith('.')):
                if other.Number == '.':
                    self.Number = other.Number
            # update type logic ***only slightly messy***
            if self.Type in ['Integer','Float']:
                if other.Type == 'Float':
                    self.Type = other.Type
                elif other.Type == 'Integer':
                    pass
                else:
                    raise Exception('Incompatible types: %s %s' % (self.Type, other.Type))
            elif self.Type in ['Character','String']:
                if other.Type == 'String':
                    self.Type = other.Type
                elif other.Type == 'Character':
                    pass
                else:
                    raise Exception('Incompatible types: %s %s' % (self.Type, other.Type))
            # description is kept
            return True
        else:
            return False

    def __lt__(self, other):
        if self.ID == 'GT':
            return True
        elif other.ID == 'GT':
            return False
        else:
            return self.ID < other.ID

    def __repr__(self):
        return '##FORMAT=<ID=%s,Number=%s,Type=%s,Description="%s">' % (self.ID,self.Number,self.Type,self.Description)


# class for variant merging
class Variant(object):
    def __init__(self,data,samplename):
        self.CHROM = data[0]
        self.POS = data[1]
        self.ID = data[2]
        self.REF = data[3]
        self.ALT = data[4]
        self.QUAL = data[5]
        self.FILTER = '.'
        self.INFO = '.'
        self.samples = { samplename: dict(zip(data[8].split(':'),data[9].split(':'))) }
        self.support = 1
        self.qualfrom = samplename
        return

    def __lt__(self,other):
        return (self.CHROM, int(self.POS)) < (other.CHROM, int(other.POS))

    def location(self):
        return self.CHROM+':'+self.POS

    def extend(self,other):
        if self.REF == other.REF and self.ALT == other.ALT and \
            (self.ID == other.ID or self.ID.startswith('rs') != other.ID.startswith('rs')):
            # update ID
            if other.ID.startswith('rs'):
                self.ID = other.ID
            # check sample names don't overlap
            try:
                assert set(self.samples.keys()).isdisjoint(set(other.samples.keys()))
            except AssertionError:
                print self.samples
                print other.samples
                raise Exception('same samples in files to be merged?')
            # merge
            if float(other.QUAL) > float(self.QUAL):
                self.QUAL = other.QUAL
                self.qualfrom = other.samples.keys()[0]
            self.samples.update(other.samples)
            self.support += 1  # increment supporting evidence
            return True
        else:
            return False

    def __repr__(self):
        return '\t'.join([self.CHROM, self.POS, self.ID, self.REF, self.ALT, \
            self.QUAL, self.FILTER, 'QUALFROM='+self.qualfrom])

    def printline(self,formatorder,sampleorder):
        formatvalues = []
        for s in sampleorder:
            formats = []
            for f in formatorder:
                try:
                    formats.append(self.samples[s][f])
                except KeyError:
                    formats.append('{nope}'.format(nope="./." if f=='GT' else '.'))
            formatvalues.append(':'.join(formats))
        return '\t'.join([ repr(self),":".join(formatorder)] + formatvalues)


# class dor header merging
class Header(object):
    def __init__(self,fi):
        self.format = {}
        self.version = None
        with gzip.open(fi) as fh:
            for line in fh:
                if line.startswith('##FORMAT'):  ## merge
                    ft = Format(line)
                    self.format[ft.ID] = ft
                elif line.startswith('##fileformat'):
                    self.version = line.rstrip()
                elif not line.startswith('#'):
                    break
        return

#    def __repr__(self):
#        #return self.format
#        return "####FAKEHEADER###"

    def extend(self, other):
        for f in other.formats.keys():
            try:
                assert self.formats[f].update(other.formats[f])
            except KeyError:
                self.formats[f] = other.formats[f]
            except AssertionError:
                raise
        self.version = sorted([self.version,other.version])[-1]

    def formatorder(self):
        return ['GT']
        #return [ x.ID for x in sorted(self.format.values()) ]


# reads slices from BED file
def readSlices(fi,slen=1000000):
    '''returns slices of given size from BED or DICT file'''
    with open(fi) as fh:
        for line in fh:
            # DICT/BED
            if line.startswith('@'):
                m = re.match('@SQ\s+SN:(\w+)\s+LN:(\d+)',line)
                if m:
                    f = [ m.group(1), 0, m.group(2) ]
                else:
                    continue
            else:
                f = line.split()
            # make segments
            for r in range(int(f[1]),int(f[2]),slen):
                yield (f[0],r,min([r+slen,int(f[2])]))


# MAIN
if __name__ == "__main__":
    from optparse import OptionParser
    usage = "USAGE: mergevcf.py -b <BED> -o <outfile.vcf> <sample1.vcf.gz> <sample2.vcf.gz>..."
    parser = OptionParser(version="mergevcf.py v1.0", usage=usage)
    #   configuration files
    parser.add_option("-o", dest="outfile",metavar="STRING", help="Merged .VCF output")
    parser.add_option("-s", dest="bedfile",metavar="STRING", help="DICT/BED file (segments to merge)")
    parser.add_option("-e", dest="evidence",metavar="INT", default=2, type=int, help="minimal evidence [2]")
    parser.add_option("-f", dest="filter",metavar="STRING,STRING,...", default='LowQual', help="exclusion filter [LowQual]")
    (options, args) = parser.parse_args()
    if len(args) < 2 or not options.outfile or not options.bedfile:
        print >> sys.stderr, __doc__
        parser.print_help()
        sys.exit(1)

    # open indexed vcfs and read headers
    vcfs = {}
    headers = None
    sampleorder = []
    for i, v in enumerate(args):
        samplename = None
        m = re.search('(?:\w{8}\.)?(..)\.',v[v.rfind('/')+1:])
        try:
            assert m
        except:
            raise
            samplename = 'sample'+str(i)
        else:
            samplename = m.group(1)
        # filehandles
        try:
            vcfs[samplename] = tabix.open(v)
        except tabix.TabixError:
            raise Exception('Cannnot open %s Is file present and tabix indexed?' % v)
        sampleorder.append(samplename)
        # headers
        try:
            headers.extend(Header(v))
        except AttributeError:
            headers = Header(v)
        except:
            raise

    # print header
    with open(options.outfile,'w') as ofh:
        print >> ofh, headers.version
        print >> ofh, '##INFO=<ID=QUALFROM,Number=1,Type=String,Description="Source of quality score (max)">'
        print >> ofh, '\n'.join(map(repr,sorted(headers.format.values())))
        print >> ofh, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + '\t'.join(sampleorder)


        stats = Counter()
        # merge and print for each slice
        for s in readSlices(options.bedfile):
            variants = {}
            sys.stderr.write('\rMerging slice '+s[0]+':'+str(s[1])+'-'+str(s[2]))
            for sample in sampleorder:
                slicevar = vcfs[sample].query(s[0], int(s[1]), int(s[2]))
                for var in slicevar:
                    # filter LowQual
                    if var[6] in options.filter.split(','):
                        continue
                    stats[sample] += 1
                    v = Variant(var,sample)
                    try:
                        success = variants[v.location()].extend(v)
                        assert success
                    except KeyError:
                        variants[v.location()] = v  # didnt exist
                    except AssertionError:
                        print >> sys.stderr, '\nWARNING incompatible variants (keeping higher scoring):'
                        ## keeping variant with higher quality score
                        stats['_CONFLICTING'] += 1
                        if float(v.QUAL) > float(variants[v.location()]):
                            variants[v.location()] = v
                            print >> sys.stderr, '\t', v
                            print >> sys.stderr, "DISCARDED:"
                            print >> sys.stderr, repr(variants[v.location()])
                        else:
                            print >> sys.stderr, repr(variants[v.location()])
                            print >> sys.stderr, "DISCARDED:"
                            print >> sys.stderr, '\t', v


            # print variants in slice
            for k, mvar in sorted(variants.items()):
                stats['_GRANDTOTAL'] += 1
                if mvar.support >= options.evidence:
                    stats['_PASSED'] += 1
                    print >> ofh, mvar.printline(headers.formatorder(),sampleorder)
        # stats
        sys.stderr.write('\r')
        for k in sorted(stats.keys()):
            sys.stderr.write('  '+k+':'+str(stats[k]))
        sys.stderr.write('\n')


