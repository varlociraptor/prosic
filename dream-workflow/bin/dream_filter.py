#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
import os
import sys
import vcf

__author__ = "Alexander Schoenhuth"

usage = """%prog <dream-vcf-file> <filtered-dream-vcf-file> <minprob> <minvaf> <pass>

	<dream-vcf-file> 	   dream-compatible vcf
        <filtered-dream-vcf-file>  dream-compatible vcf, filtered by probability and VAF
        <minprob>                  minimum probability
        <minvaf>                   minimum variant allele frequency (VAF)
        <passtrue>                     if set to 1, picks only 'PASS' variants

        The format of <dream-vcf-file> is mandatory for the dream challenge. 

	1) This VCF will ONLY contain the SOMATIC calls, and, 
	2) will follow the VCF 4.1 format. 

Note that the challenge provides a program to test the correctness of the 
format. 
"""

def main():
    
    parser = OptionParser(usage=usage)
    (options, args) = parser.parse_args()

    if (len(args)!=5):
        parser.print_help()
        return 1

    dream_vcf_filename = os.path.abspath(args[0])
    filtered_dream_vcf_filename = os.path.abspath(args[1])
    minprob = float(args[2])
    minvaf = float(args[3])
    passtrue = (args[4] == '1')

#    print(dream_vcf_filename, filtered_dream_filename)
    dream_records = open(dream_vcf_filename, 'r').readlines()
    filtered_records = open(filtered_dream_vcf_filename, 'w')

    for recstr in dream_records:
        
        rec = recstr.strip().split()

        if rec[0][0] == '#':
            filtered_records.write(recstr)
            continue

        prob = float(rec[7].strip().split(';')[4].split('=')[1])
        vaf =  float(rec[7].strip().split(';')[3].split('=')[1])
        
        if passtrue:
            if rec[6] != 'PASS':
                continue

        #print(prob, vaf, rec[6])
        if prob < minprob or vaf < minvaf:
            continue

        filtered_records.write(recstr)

if __name__ == '__main__':
    sys.exit(main())
