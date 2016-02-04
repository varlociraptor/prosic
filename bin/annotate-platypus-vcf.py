#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
import os
import sys
import vcf

__author__ = "Alexander_Schoenhuth"

usage = """%prog [options] <platypus-vcf-file> <new-vcf-file> 

	<vcf-file> 	tabix-indexed VCF file 
	<new-vcf-file> 	annoted VCF file

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
        if gt_tumour != '0/0' and gt_control == '0/0':
                return 'SOMATIC', 1.0
        if gt_control != '0/0':
                return 'GERMLINE', 1.0
        return 'ABSENT', 1.0

def determinePlatypusVAFs (gt_tumour, gt_control):
        
        if gt_tumour == '0/0':
                c_vaf = 0.0
        elif gt_tumour in ['0/1','1/0']:
                c_vaf = 0.5
        elif gt_tumour == '1/1':
                c_vaf = 1.0
        else:
                c_vaf = 'UNKNOWN'
        if gt_control == '0/0':
                h_vaf = 0.0
        elif gt_control in ['0/1','1/0']:
                h_vaf = 0.5
        elif gt_control == '1/1':
                h_vaf = 1.0
        else:
                h_vaf = 'UNKNOWN'

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

        for vcf_record in vcf_reader:

                if vcf_record.FILTER != 'PASS':
                        continue
                if len(vcf_record.ALT) > 1:
                        continue
                
                cancer_gt, healthy_gt = vcf_record.samples[0]['GT'], vcf_record.samples[1]['GT']

                if cancer_gt == './.' or healthy_gt == './.':
                        continue
                
                call, post_prob_call = makeDKFZPlatypusCall(cancer_gt, healthy_gt)
                c_vaf, h_vaf = determinePlatypusVCFs(cancer_gt, healthy_gt)

#		vcf_record = vcf_reader.next() 

		vcf_record.INFO['CALL'] 		= call 
		vcf_record.INFO['POSTERIOR_PROB'] 	= post_prob_call
		vcf_record.INFO['MAP_HEALTHY_VAF'] 	= h_vaf
		vcf_record.INFO['MAP_CANCER_VAF'] 	= c_vaf
                
		vcf_writer.write_record(vcf_record)

#	for vcf_record in vcf_reader:
		vcf_writer.write_record(vcf_record)
			
			 

if __name__ == '__main__':
	sys.exit(main())
