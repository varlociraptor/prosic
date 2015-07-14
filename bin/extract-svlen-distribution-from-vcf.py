#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
import os
import vcf
import sys
import operator
from collections import defaultdict

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__))[:-3] + 'python')
from Indel import *

__author__ = "Louis Dijkstra"

usage = """%prog [options] <vcf-file> 

	<vcf-file> 		original VCF file

Outputs the data for the length distribution of the indels in the given VCF file. 
The output is organized in two column (tab-seperated): 

	x_1	c_1
	x_2	c_2
	...	...
	x_n	c_n

where x_1 is the minimal value found and x_n is the maximum value found. (Note: x_{i+1} = x_i + 1). 
c_i is the count for x_i. 
"""

class CountData:

	def __init__(self):
		self.n = 0 # number of data points added
		self.min_value = float('Inf')
		self.max_value = float('-Inf')
		self.count = defaultdict(int)	

	def add(self, value, count):
		self.n += 1 
		if value < self.min_value: 
			self.min_value = value 
		if value > self.max_value:
			self.max_value = value 
		self.count[value] += count

	def print(self, output_stream=sys.stdout):
		if self.n == 0:
			return 0
		for i in range(self.min_value, self.max_value + 1):
			print("%d\t%d"%(i, self.count[i]), file=output_stream) 
		

def main():

	parser = OptionParser(usage=usage)
	parser.add_option("-o", action="store", dest="outputfilename", default=None,  
				help="Data is stored under this name.")
	parser.add_option("-v", action="store_true", dest="verbose", default=False,
                      		help="Verbose. Prints regularly how many variants have been processed.")
	(options, args) = parser.parse_args()
	
	if (len(args)!=1):
		parser.print_help()
		return 1

	vcf_filename 		= os.path.abspath(args[0])
	vcf_reader 		= vcf.Reader(open(vcf_filename))
	
	svlen_distr = CountData()

	# Walk through all records in the VCF file
	n = 0
	for vcf_record in vcf_reader: 
		n += 1 
		if options.verbose and n % 10000 == 0: 
			print("Processed %d variants"%n) 
		is_deletion = isDeletion(vcf_record)
		is_insertion = isInsertion(vcf_record) 
		if not is_deletion and not is_insertion: # no indel? 
			continue 	

		length = returnIndelLength(vcf_record)
		if is_deletion: 
			svlen_distr.add(-1 * length, 1) 
		else: 
			svlen_distr.add(length, 1) 
		
	if options.outputfilename != None: 
		svlen_distr.print(output_stream = open(options.outputfilename, 'w'))
	else: 
		svlen_distr.print()

if __name__ == '__main__':
	sys.exit(main())

