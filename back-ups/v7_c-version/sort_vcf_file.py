#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
import os
import vcf
import sys
import operator

__author__ = "Louis Dijkstra"

usage = """%prog [options] <vcf-file> <output-vcf-file>

Sorts all records in a VCF file by chromosome and then by position. 
 
	<vcf-file> 		original VCF file
	<output-vcf-file>	name where the sorted VCF records will be stored	
"""

def returnChromosome (chromosome):
	"""Returns an integer denoting the chromosome."""
	if len(chromosome) > 3 and chromosome[:3] == 'chr':
		chromosome = chromosome[3:]
	if chromosome == 'X' or chromosome == 'x':
		return 23
	if chromosome == 'Y' or chromosome == 'y':
		return 24	
	return int(chromosome)

def main():

	parser = OptionParser(usage=usage)
	parser.add_option("-x", action="store", dest="chromosome", default=None, type=int, 
						help="Processes only this chromosome. (Default = all are processed)")
	(options, args) = parser.parse_args()
	
	if (len(args)!=2):
		parser.print_help()
		return 1

	# list of chromosomes to be processed 
	list_of_chromosomes = range(1,25) # default: all chromosomes
	if options.chromosome is not None:
		list_of_chromosomes = [options.chromosome]
	
	vcf_filename 		= os.path.abspath(args[0])
	vcf_output_filename	= os.path.abspath(args[1])

	vcf_reader 		= vcf.Reader(open(vcf_filename))
	vcf_writer 		= vcf.Writer(open(vcf_output_filename, 'w'), vcf_reader)

	list_vcf_records = [[] for i in range(24)] # one list for every chromosome

	# Walk through all records in the VCF file
	for vcf_record in vcf_reader: 
		chromosome = returnChromosome(vcf_record.CHROM)
		if not chromosome in list_of_chromosomes: # not the chromosome we are currently interested in.		
			continue
		list_vcf_records[chromosome - 1].append(vcf_record)
		
	
	# sort the chromosomes by position
	for chromosome in list_of_chromosomes: 
		list_vcf_records[chromosome - 1].sort(key=operator.attrgetter('POS')) 

	# write them to the new vcf file:
	for chromosome in list_of_chromosomes: 
		for vcf_record in list_vcf_records[chromosome - 1]:
			vcf_writer.write_record(vcf_record)
		
	vcf_writer.close()

if __name__ == '__main__':
	sys.exit(main())

