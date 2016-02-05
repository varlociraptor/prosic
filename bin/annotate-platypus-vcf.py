#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
import os
import sys
import vcf

__author__ = "Alexander_Schoenhuth"

usage = """%prog [options] <platypus-vcf-file> <new-vcf-file> 

	<vcf-file> 	tabix-indexed VCF file 
	<new-vcf-file> 	annotated VCF file

Annotates the given VCF file with the results generated by the DKFZ
calling routine for platypus.
"""

def makeCall (p_somatic, p_germline, p_not_present):
	"""Makes the call on the basis of the posterior probabilities."""
	if p_somatic > p_germline and p_somatic > p_not_present:
		return 'SOMATIC', p_somatic
	if p_germline > p_somatic and p_germline > p_not_present:
		return 'GERMLINE', p_germline
	return 'ABSENT', p_not_present	

def makePlatypusCall (gt_tumour, gt_control):
	"""Makes the call on the basis of the tumour-control genotypes."""
        if gt_tumour in ['0/1','1/0','1/1'] and gt_control == '0/0':
                return 'SOMATIC', 1.0
        elif gt_control in ['0/1','1/0','1/1'] and gt_tumour:
                return 'GERMLINE', 1.0
        elif gt_tumour == '0/0' and gt_control == '0/0':
                return 'ABSENT', 1.0
        else:
#                print(gt_tumour, gt_control, file=sys.stderr)
                return 'UNKNOWN', 'NONE'
        
def determinePlatypusVAFs (gt_tumour, gt_control):
        
        if gt_tumour == '0/0':
                c_vaf = 0.0
        elif gt_tumour in ['0/1','1/0']:
                c_vaf = 0.5
        elif gt_tumour == '1/1':
                c_vaf = 1.0
        else:
                c_vaf = 'NONE'
        if gt_control == '0/0':
                h_vaf = 0.0
        elif gt_control in ['0/1','1/0']:
                h_vaf = 0.5
        elif gt_control == '1/1':
                h_vaf = 1.0
        else:
                h_vaf = 'NONE'

        return c_vaf, h_vaf
                        
                        
def main():

	parser = OptionParser(usage=usage)
	parser.add_option("-v", action="store_true", dest="verbose", default=False,
                      		help="Verbose. Prints regularly how many variants have been processed.")
	(options, args) = parser.parse_args()
	
	if (len(args)!=2):
		parser.print_help()
		return 1
	
	vcf_filename 		= os.path.abspath(args[0])
	new_vcf_filename 	= os.path.abspath(args[1])

	vcf_reader 		= vcf.Reader(open(vcf_filename))
	vcf_writer 		= vcf.Writer(open(new_vcf_filename, 'w'), vcf_reader)

        i = 0
        for vcf_record in vcf_reader:

                #print(vcf_record.CHROM, vcf_record.POS, vcf_record.FILTER, file=sys.stderr)
                if len(vcf_record.FILTER) > 0:
                        #print("didn't PASS", file=sys.stderr)
                        continue
                i += 1
                if i % 1000 == 0:
                        print(i, "good variants so far", file=sys.stderr)
                if len(vcf_record.ALT) > 1:
#                        print("multilalelic site", file=sys.stderr)
                        continue
                
                cancer_gt, healthy_gt = vcf_record.samples[0]['GT'], vcf_record.samples[1]['GT']
                
                if not (cancer_gt and healthy_gt):
                        print("One genotype = ./.:", vcf_record.CHROM, vcf_record.POS, file=sys.stderr)
                        continue
                
                call, post_prob_call = makePlatypusCall(cancer_gt, healthy_gt)
                if call == "UNKNOWN":
                        print("UNKNOWN: ", vcf_record.CHROM, vcf_record.POS, cancer_gt, healthy_gt, file=sys.stderr)
                        
                c_vaf, h_vaf = determinePlatypusVAFs(cancer_gt, healthy_gt)

		vcf_record.INFO['CALL'] 		= call 
		vcf_record.INFO['POSTERIOR_PROB'] 	= post_prob_call
		vcf_record.INFO['MAP_HEALTHY_VAF'] 	= h_vaf
		vcf_record.INFO['MAP_CANCER_VAF'] 	= c_vaf
                
		vcf_writer.write_record(vcf_record)

#	for vcf_record in vcf_reader:
		vcf_writer.write_record(vcf_record)
			
        print(i, "good variants in total", file=sys.stderr)

if __name__ == '__main__':
	sys.exit(main())
