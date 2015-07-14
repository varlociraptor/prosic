#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
import os
import vcf
import sys

__author__ = "Louis Dijkstra"

from Variant.Indel import *

usage = """%prog [options] <vcf-file> <output-vcf-file>

Extracts all indels with just one alternative from the VCF file <vcf-file> and 
stores them into a new VCF file <output-vcf-file>. 
Indels on the sex chromosomes are discarded.  
 
	<vcf-file> 		original VCF file
	<output-vcf-file>	name of VCF file where the deletions will be stored.	
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

def main():

	parser = OptionParser(usage=usage)
	parser.add_option("--deletions-only", action="store_true", dest="deletions_only", default=False, 
				  		help="Only deletions are stored and insertions are discarded.")
	parser.add_option("--insertions-only", action="store_true", dest="insertions_only", default=False, 
				  		help="Only insertions are stored and deletions are discarded.")
	parser.add_option("-m", action="store", dest="min_length", default=0, type=int,
				  		help="Minimal length of a deletion to be considered. (Default = 0 bp; no threshold)")
	parser.add_option("-x", action="store", dest="autosome", default=None, type=int, 
						help="Processes only this autosome. (Default = all are processed)")
	(options, args) = parser.parse_args()
	
	if (len(args)!=2):
		parser.print_help()
		return 1

	# list of autosomes to be processed 
	list_of_autosomes = range(1,23) # default: all autosomes
	if options.autosome is not None:
		list_of_autosomes = [options.autosome]
	
	vcf_filename 		= os.path.abspath(args[0])
	vcf_output_filename	= os.path.abspath(args[1])

	vcf_reader 		= vcf.Reader(open(vcf_filename))
	vcf_writer 		= vcf.Writer(open(vcf_output_filename, 'w'), vcf_reader)

	# Walk through all records in the VCF file
	for vcf_record in vcf_reader:  
		if not returnAutosome(vcf_record.CHROM) in list_of_autosomes: # not the autosome we are currently interested in.		
			continue

		if len(vcf_record.ALT) != 1: # records with several alternatives are ignored.
			continue 
		
		is_deletion 	= isDeletion(vcf_record) 
		is_insertion 	= isInsertion(vcf_record)
		if not (is_deletion or is_insertion): # not an indel
			continue
		if options.only_deletions and is_insertion:
			continue
		if options.only_insertions and is_deletion:
			continue
		if returnIndelLength(vcf_record) < options.min_length: # below the minimal length of an indel.
			continue 
		
		vcf_writer.write_record(vcf_record)
		
	vcf_writer.close()


if __name__ == '__main__':
	sys.exit(main())
