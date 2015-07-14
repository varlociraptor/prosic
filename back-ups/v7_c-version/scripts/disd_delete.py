#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function, division
from optparse import OptionParser
import sys
import os
from pysam import Samfile
import math
import scipy.stats

__author__ = "Alexander Schoenhuth"

usage = """%prog [options]

Reads BAM file and returns median and and median absolute deviation of the insert size distribution
of its reads as robust estimators for mean and variance of the normal distribution which reflects
the insert size distribution. Mean is calculated as median whereas standard deviation is calculated
as 1.4826 * median absolute deviation

"""

def isUnique(read_aln):
    
    """return True if read has probability one (== YP:f:1), False otherwise"""
    
    ysindex = None

    for tag in read_aln.tags:
        if tag[0] == 'YP':
            ypindex = read_aln.tags.index(tag)

    if ypindex == None:
        return False
    else:
        return read_aln.tags[ypindex][1] == 1 
    
    #ysindex = [read_aln.tags.index(tag) for tag in read_aln.tags if tag[0] == 'YS'][0]

def computeMedian(sizedict):
    
    overallcount = 0
    intlengths = sizedict.keys()
    intlengths.sort()
    for l in intlengths:
        overallcount += sizedict[l]

    accucount = 0
    for l in intlengths:
        accucount += sizedict[l]
        if accucount > float(overallcount)/2:
            median = l
            break

    return median

def main():
    
    parser = OptionParser(usage=usage)
    parser.add_option("-s", action="store_true", dest="sam_input", default=False,
                      help="Input is in SAM format instead of BAM format")
    parser.add_option("-S", action="store_true", dest="sam_output", default=False,
                      help="Write output in SAM format instead of in BAM format")
    parser.add_option("-l", action="store", type="int", dest="read_length", default=100,
                      help="read length, default is 100")
    (options, args) = parser.parse_args()
    if (len(args)!=0) or (os.isatty(0)):
        parser.print_help()
        sys.exit(1)

    inputfile = Samfile('-', 'r' if options.sam_input else 'rb')
#    outfile = Samfile('-', 'wh' if options.sam_output else 'wb', template=inputfile)
    
    sizedict = {}

    counter = 0
    for read_aln in inputfile:

        counter += 1
        if counter % 100000 == 0:
            print("having processed %d read alignments" % (counter), file=sys.stderr)

#        if not isUnique(read_aln):
#            continue
#        if not read_aln.is_read2:
#            continue # we only want to count insert size once for each read
        
        intsize = read_aln.isize - 2*options.read_length
        if intsize == 0 or intsize < 0: # < 0 iff read_aln.is_read1 so no double counting
            continue
        if sizedict.has_key(intsize):
            sizedict[intsize] += 1
        else:
            sizedict[intsize] = 1
    
#    print(sizedict)
    median = computeMedian(sizedict)
    
    intlengths = sizedict.keys()
    intlengths.sort()

    deviationdict = {}
    for l in intlengths:
        if deviationdict.has_key(abs(l-median)):
            deviationdict[abs(l-median)] += sizedict[l]
        else:
            deviationdict[abs(l-median)] = sizedict[l]
    
    medianabsdeviation = computeMedian(deviationdict)
    outstring = "Median: %d, median absolute deviation: %d" % (median, medianabsdeviation)
    print(outstring)
    
    sys.exit(1)

if __name__ == '__main__':
    sys.exit(main())
