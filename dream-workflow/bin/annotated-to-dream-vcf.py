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

    i, j, k = 0, 0, 0
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
        #if len(filterr) > 0: # FILTER other than ['PASS']
        #    continue
        if filterr == None or len(filterr) == 0:
            fstr = 'PASS'
        else:
            fstr = ''
            for x in filterr[:-1]:
                fstr += "%s;" % (x)
            fstr += filterr[-1]
            
        info = vcf_record.INFO

        if not info['SVTYPE'] in ['INS','DEL'] or str(alt)[0] == '<':
            #print(info['SVTYPE'],file=sys.stderr)
            #print("continuing, found structvar...",file=sys.stderr)
            j += 1
            continue

        #print(str(alt), file=sys.stderr)
        
        svlen = len(alt) - len(ref)
        endpos = pos + len(ref) - 1

        #print(chrom, pos, idd, ref, alt, qual, filterr, info['CALL'], svlen, endpos)
        #print(info, file=sys.stderr)
        if not 'CALL' in info or not 'POSTERIOR_PROB' in info:
            k += 1
            #print("continuing, no call specification, no probability...",file=sys.stderr)
            continue
        #print(info['CALL'],file=sys.stderr)
        somatic = (info['CALL'][0] == 'SOMATIC')
        
        if not somatic:
            continue
        else:
            #print("Somatic!!!", file=sys.stderr)
            vaf = float(info['MAP_CANCER_VAF'][0])
            prob = float(info['POSTERIOR_PROB'][0])
            #if prob < minmapp:
            #    continue

        # writing original genotypes, if available

        hgt = vcf_record.genotype('normal')['GT']
        if not hgt == '0/0':
            continue
#        if tgt in ['1/0','0/1']:
#            tgt = 0.5
#        elif tgt == '0/0':
#            tgt = 0.0
#        elif tgt == '1/1':
#            tgt = 1.0
#        else:
#            continue
#        tad = vcf_record.genotype('tumor')['AD']
#        
#        hgt = vcf_record.genotype('normal')['GT']
#        if hgt in ['1/0','0/1']:
#            hgt = 0.5
#        elif hgt == '0/0':
#            hgt = 0.0
#        elif hgt == '1/1':
#            hgt = 1.0
#        else:
#            continue
#        had = vcf_record.genotype('normal')['AD']

        i += 1
        if (i % 100 == 0):
            print("Found %d somatic calls so far" % (i)) 

#        dream_vcf.write('%s\t%d\t%s\t%s\t%s\t.\t%s\tSOMATIC;END=%d;SVLEN=%d;VAF=%.6f;PROB=%.8f;TGT=%.2f;TAD=%d,%d;HGT=%.2f;HAD=%d,%d\n' % (chrom, pos, idd, ref, alt, fstr, endpos, svlen, vaf, prob, tgt, int(tad[0]), int(had[1]), hgt, int(had[0]), int(had[1])))
        dream_vcf.write('%s\t%d\t%s\t%s\t%s\t.\t%s\tSOMATIC;END=%d;SVLEN=%d;VAF=%.6f;PROB=%.8f\n' % (chrom, pos, idd, ref, alt, fstr, endpos, svlen, vaf, prob))

    print("Found %d structural variants" % (j), file=sys.stderr)
    print("%d indels not annotated" % (k), file=sys.stderr)
    print("Overall %d somatic calls" % (i), file=sys.stderr)
        
if __name__ == '__main__':
    sys.exit(main())

    
