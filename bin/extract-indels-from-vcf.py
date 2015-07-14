#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
import os
import vcf
import sys

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__))[:-3] + 'python')
from Indel import *

__author__ = "Louis Dijkstra"

usage = """%prog [options] <vcf-file> <output-vcf-file>

	<vcf-file> 		original VCF file
	<output-vcf-file>	name of VCF file where the deletions will be stored.	

Extracts all indels from the VCF file <vcf-file> and 
stores them into a new VCF file <output-vcf-file>. 
"""

def main():

	parser = OptionParser(usage=usage)
	parser.add_option("-v", action="store_true", dest="verbose", default=False,
                      		help="Verbose. Prints regularly how many variants have been processed.")
	(options, args) = parser.parse_args()
	
	if (len(args)!=2):
		parser.print_help()
		return 1

	vcf_filename 		= os.path.abspath(args[0])
	vcf_output_filename	= os.path.abspath(args[1])

	vcf_reader 		= vcf.Reader(open(vcf_filename))
	vcf_writer 		= vcf.Writer(open(vcf_output_filename, 'w'), vcf_reader)

	# Walk through all records in the VCF file
	total, n_indels = 0, 0 
	for vcf_record in vcf_reader:  
		total += 1 

		if isDeletion(vcf_record) or isInsertion(vcf_record): 
			vcf_writer.write_record(vcf_record)
			n_indels += 1 

		if options.verbose and total % 1000 == 0: 
			print("Processed %d variants of which %d were indels"%(total, n_indels))
		
	vcf_writer.close()


if __name__ == '__main__':
	sys.exit(main())
