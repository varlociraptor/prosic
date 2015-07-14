#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
import os
import sys
import math

from summarize_results import *

__author__ = "Louis Dijkstra"

usage = """%prog <results-filename>

Plots the cancer VAF MLEs against the true cancer VAFs. 

	<results-filename> 	File with the somatic mutation results. 
"""



def main():

	parser = OptionParser(usage=usage)
	parser.add_option("-a", action="store", dest="length_start", default=None, type=int,
						help="Minimal length of the indels taken into account. (Default = unbounded)")
	parser.add_option("-b", action="store", dest="length_end", default=None, type=int,
						help="Maximal length of the indels taken into account. (Default = unbounded)")
	(options, args) = parser.parse_args()
	
	if (len(args)!=1):
		parser.print_help()
		return 1

	results_filename = os.path.abspath(args[0])
	results_file = open(results_filename, 'r')

	results = [] # list with all the results
	c_vafs = [] # a list of all the different cancer VAFs present in the file

	# walk through all the results in the file
	for line in results_file: 
		result = ResultsInstance(line)
		if result.chromosome == 'X' or result.chromosome == 'X ' or result.chromosome == 'Y' or result.chromosome == 'Y ' or result.chromosome == '': # sex chromosomes are not taken into account
			continue
		
		if result.global_max_exists == 0:
			continue 

		# check length range
		if options.length_start != None:
			if result.length < options.length_start:
				continue
		if options.length_end != None:
			if result.length > options.length_end:
				continue

		if result.true_c_vaf not in c_vafs: 
			c_vafs.append(result.true_c_vaf)
		results.append(result) 

	c_vafs.sort()

	n_different_c_vafs = len(c_vafs)

	mles = []
	for i in range(n_different_c_vafs):
		mles.append([])

	for result in results:
		mles[c_vafs.index(result.true_c_vaf)].append(result.returnCancerVAF_MLE())

	for i in range(len(mles)):
		mles[i].sort()

	print mles
		

	


				
	
if __name__ == '__main__':
	sys.exit(main())
