#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
import os
import vcf
import sys
import operator

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__))[:-3] + 'python')
from Indel import *

__author__ = "Louis Dijkstra"

usage = """%prog [options] <vcf-file> <output-vcf-file>
 
	<vcf-file> 			original VCF file
	<output-vcf-file>	filtered VCF records will be stored here
	
Filters a given VCF file. Type '%prog -h' for the filter options. 	
"""


def determineType (ref, alt): 
	"""Determine whether the VCF record is a SNP, deletion or insertion
	   on the basis of the reference and alternative sequence.

	   Returns 3 booleans: whether its a 1) SNP, 2) deletion, or 3) insertion."""
	len_ref = len(ref)
	len_alt = len(alt)
	
	if len_ref == 1 and len_alt == 1:
		return True, False, False
	if len_ref > 1 and len_alt == 1:
		return False, True, False			
	if len_ref == 1 and len_alt > 1:
		return False, False, True		
	return False, False, False			
				

def main():

	parser = OptionParser(usage=usage)
	parser.add_option("--deletions", action="store_true", dest="deletions", default=False, 
				help = "Filter for deletions")
	parser.add_option("--NOT", action="store_true", dest="negation", default=False, 
				help = "Outputs those records that do NOT satisfy the given constraints. E.g., '--deletions --NOT' returns all records that are not deletions.") 
	parser.add_option("--insertions", action="store_true", dest="insertions", default=False, 
				help = "Filter for insertions")
	parser.add_option("--snps", action="store_true", dest="snps", default=False, 
				help = "Filter for SNPs")
	parser.add_option("--use-ref-alt", action="store_true", dest="use_ref_alt", default=False, 
				help = "Uses the REF and ALT field to determine whether an entry is an indel or not. Otherwise the INFO field is used.")
	parser.add_option("-m", action="store", dest="min_length", default=0, type=int, 
				help="Minimal length of the indels (Default = 0; i.e., all indels are outputed)")
	parser.add_option("-v", action="store_true", dest="verbose", default=False,
				help="Verbose. Prints regularly how many variants have been processed.")
	parser.add_option("-x", action="store", dest="chromosome", default=None, 
				help="Processes only this chromosome. (Default = all are processed)")
	(options, args) = parser.parse_args()
	
	if (len(args)!=2):
		parser.print_help()
		return 1
	
	vcf_filename 		= os.path.abspath(args[0])
	vcf_output_filename	= os.path.abspath(args[1])

	vcf_reader 		= vcf.Reader(open(vcf_filename))
	vcf_writer 		= vcf.Writer(open(vcf_output_filename, 'w'), vcf_reader)

	n = 0 
	for vcf_record in vcf_reader:
		
		n += 1
		if options.verbose and n % 10000 == 0: 
			print('Processed %d variants'%n)
		
		satisfies_filter = True
		if options.chromosome != None and options.chromosome != vcf_record.CHROM: 
			satisfies_filter = False

		if satisfies_filter:
			satisfies_filter 	= False

			is_deletion 		= isDeletion(vcf_record) 
			is_insertion 		= isInsertion(vcf_record)
			is_snp 			= isSNP(vcf_record) 

			if options.use_ref_alt: # use the REF and ALT field to classify
				is_snp, is_deletion, is_insertion = determineType (vcf_record.REF, vcf_record.ALT[0])

			if options.snps and is_snp: 
				satisfies_filter = True
			if options.deletions and is_deletion: 
				if returnIndelLength(vcf_record) >= options.min_length:
					satisfies_filter = True
			if options.insertions and is_insertion: 
				if returnIndelLength(vcf_record) >= options.min_length:
					satisfies_filter = True

		if not options.negation and satisfies_filter:
			vcf_writer.write_record(vcf_record)
		if options.negation and not satisfies_filter:
			vcf_writer.write_record(vcf_record)	
		
	vcf_writer.close()

if __name__ == '__main__':
	sys.exit(main())

