#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
import os
import sys
import vcf

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__))[:-3] + 'python')

from Indel import * 

__author__ = "Louis Dijkstra"

usage = """%prog <true-vcf-file> 

	<true-vcf-file> 	tabix-indexed VCF file 

Extracts the true calls (somatic/germline/absent) from the VCF file
containing the ground truth. 
"""

def returnAutosome (autosome):
	"""Returns an integer denoting the chromosome."""
	if len(autosome) > 3 and autosome[:3] == 'chr':
		autosome = autosome[3:]
	if autosome == 'X' or autosome == 'x':
		return 23
	if autosome == 'Y' or autosome == 'y':
		return 24	
	return int(autosome)

def determineVAF(gt_nums):
	if gt_nums == '1|1':
		return 1.0 
	elif gt_nums == None:
		return 0.0 
	return 0.5
			
def returnTruth(true_vcf_record):
	"""Returns the true VAFs of the healthy and cancer cells and the classification (somatic/germline)"""
	# There is one control VAF and VAFs for the four cancer populations
	true_h_vaf, som1_vaf, som2_vaf, som3_vaf, som4_vaf = 0.0, 0.0, 0.0, 0.0, 0.0 
	for call in true_vcf_record.samples:
		if call.sample=='Control': 	true_h_vaf 	= determineVAF(call.gt_nums)
		if call.sample=='Som1':		som1_vaf 	= determineVAF(call.gt_nums)	
		if call.sample=='Som2':		som2_vaf 	= determineVAF(call.gt_nums)	
		if call.sample=='Som3':		som3_vaf 	= determineVAF(call.gt_nums)
		if call.sample=='Som4':		som4_vaf 	= determineVAF(call.gt_nums)
	true_c_vaf = (1/3.0)*som1_vaf + (1/3.0)*som2_vaf + (1/4.0)*som3_vaf + (1/12.0)*som4_vaf 

	true_class = 'GERMLINE'
	if true_h_vaf == 0.0:
		true_class = 'SOMATIC'

	return true_h_vaf, true_c_vaf, true_class 			

def main():

	parser = OptionParser(usage=usage)
	(options, args) = parser.parse_args()
	
	if (len(args)!=1):
		parser.print_help()
		return 1
	
	vcf_filename 		= os.path.abspath(args[0])
	vcf_reader 		= vcf.Reader(open(vcf_filename))

	for autosome in range(1,23): 
		t_deletion 		= [] 
		t_centerpoints 		= [] 
		t_length		= []
		t_h_vaf			= []
		t_c_vaf 		= [] 
		t_class 		= [] 

		try:
			vcf_record = vcf_reader.next() 

			while (returnAutosome(vcf_record.CHROM) == autosome):
				is_deletion = isDeletion(vcf_record)
				is_insertion = isInsertion(vcf_record) 
				if (not is_deletion and not is_insertion): # not an indel
					vcf_record = vcf_reader.next()  
				else: 
					length = int(returnIndelLength(vcf_record))
					t_length.append(length)
					centerpoints = [] 
					if is_deletion:	
						t_deletion.append(1) 
						if length % 2 == 0: 
							centerpoints = [vcf_record.POS + int(length / 2), vcf_record.POS + int(length / 2) + 1]
						else: 
							centerpoints = [vcf_record.POS + int(length / 2)] 
					else: 			
						t_deletion.append(0)
						centerpoints = [vcf_record.POS, vcf_record.POS + 1]
					t_centerpoints.append(centerpoints) 	
					true_h_vaf, true_c_vaf, true_class = returnTruth(vcf_record) 	
					t_h_vaf.append(true_h_vaf)
					t_c_vaf.append(true_c_vaf)
					t_class.append(true_class) 
					try:
						vcf_record = vcf_reader.next() 
					except: 
						break 
		except:
			pass 

		t_centerpoints, t_deletion, t_length, t_h_vaf, t_c_vaf, t_class = (list(t) for t in zip(*sorted(zip(t_centerpoints, t_deletion, t_length, t_h_vaf, t_c_vaf, t_class))))			
		
		for i in range(len(t_centerpoints)):
			print(autosome, '\t', end = '')
			print(t_deletion[i], '\t', end = '')
			print(t_length[i], '\t', end = '')
			print(t_h_vaf[i], '\t', end = '')
			print(t_c_vaf[i], '\t', end = '')
			print(t_class[i], '\t', end = '')
			print(t_centerpoints[i][0], end = '')
			if len(t_centerpoints[i]) > 1:
				print('\t', t_centerpoints[i][1])
			else:
				print("")
		 

if __name__ == '__main__':
	sys.exit(main())
