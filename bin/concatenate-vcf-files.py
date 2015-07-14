#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
import os
import vcf
import sys
import operator

__author__ = "Louis Dijkstra"

usage = """%prog [options] <vcf-file_1> <vcf-file_2> ... <vcf-file_m>

Concatenates the VCF files and outputs the result to the command line. 

Contatenation fails when the headers of the files differ.

NOTE: output is not sorted. See "sort-vcf-file.py".
"""

def returnHeader(vcf_file):
	"""Returns the header of the VCF file."""
	for line in vcf_file:
		if line[:6] == '#CHROM':
			return line

def main():

	parser = OptionParser(usage=usage)
	(options, args) = parser.parse_args()
	
	if (len(args) < 1):
		parser.print_help()
		return 1

	# print meta-data 
	print('##fileformat=VCFv4.1')
	print('##source=concatenate-vcf-files.py')
	
	vcf_filename = args[0]
	vcf_file = open(vcf_filename, 'r')
	header_old = returnHeader(vcf_file)
	
	print(header_old, end='')

	for line in vcf_file:
		print(line, end='')
	vcf_file.close()

	for vcf_filename in args[1:]:
		vcf_file = open(vcf_filename, 'r') 
		header_new = returnHeader(vcf_file)
		if header_old != header_new:
			print("\nERROR: headers of the VCF files differ.")
			return 1
		for line in vcf_file:
			print(line, end='')
		vcf_file.close()

if __name__ == '__main__':
	sys.exit(main())

