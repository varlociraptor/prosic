#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
import os
import sys
import vcf

__author__ = "Alexander Schoenhuth"

usage = """%prog [options] <annotated-vcf-file> <dream-vcf-file> 

	<annotated-vcf-file> 	tabix-indexed, self-annotated (by annotate-* script) VCF file 
        <dream-vcf-file> 	dream challenge compatible VCF file
        <method-name>           name as listed in the challenge
        <minmapp>               minimum a posteriori probability required for passing as somatic

        The format of <dream-vcf-file> is mandatory for the dream challenge. 

	1) This VCF will ONLY contain the SOMATIC calls, and, 
	2) will follow the VCF 4.1 format. 

Note that the challenge provides a program to test the correctness of the 
format. 
"""

def main():

    parser = OptionParser(usage=usage)
    (options, args) = parser.parse_args()

    if (len(args)!=4):
        parser.print_help()
        return 1

    annotated_vcf_filename = os.path.abspath(args[0])
    dream_vcf_filename = os.path.abspath(args[1])
    methodname = args[2]
    minmapp = float(args[3])
    
    vcf_reader = vcf.Reader(open(annotated_vcf_filename))
    dream_vcf = open(dream_vcf_filename, 'w')

    dream_vcf.write('##fileformat=VCFv4.1\n')
    dream_vcf.write('##source=%s\n' % (methodname))
    dream_vcf.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n')
    dream_vcf.write('##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Indicates if record is a somatic mutation">\n')
    dream_vcf.write('##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">\n')
    dream_vcf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')

    i = 0
    for vcf_record in vcf_reader:
        #vcf_record = vcf_reader.next()

        #print(vcf_record)
        
        chrom = vcf_record.CHROM
        pos = vcf_record.POS
        idd = vcf_record.ID
        if not idd:
            idd = '.'
        ref = vcf_record.REF
        alt = vcf_record.ALT[0]
        qual = vcf_record.QUAL
        filterr = vcf_record.FILTER
        if len(filterr) == 0:
            filterr = 'PASS'
        info = vcf_record.INFO
        
        svlen = len(alt) - len(ref)
        endpos = pos + len(ref) - 1

        #print(chrom, pos, idd, ref, alt, qual, filterr, info['CALL'], svlen, endpos)
        #print(info, file=sys.stderr)
        if not 'CALL' in info or not 'POSTERIOR_PROB' in info:
            continue
        
        somatic = (info['CALL'][0] == 'SOMATIC')
        
        if not somatic:
            continue
        else:
            posteriorprob = float(info['POSTERIOR_PROB'][0])
            if posteriorprob < minmapp:
                continue
        
        i += 1
        if (i % 100 == 0):
            print("Found %d somatic calls so far" % (i)) 

        dream_vcf.write('%s\t%d\t%s\t%s\t%s\t.\tPASS\tSOMATIC;END=%d;SVLEN=%d\n' % (chrom, pos, idd, ref, alt, endpos, svlen))

    print("Overall %d somatic calls" % (i), file=sys.stderr)
        
if __name__ == '__main__':
    sys.exit(main())

    
