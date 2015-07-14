#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
import os
import sys
import numpy as np
from itertools import islice

__author__ = "Louis Dijkstra"

usage = """%prog <observations-file> 

	<observations-file> 	observations file generated with
				'extract-observations.py' 

Returns robust estimates of the mean and standard deviation of the insert size
distribution in both the healthy and the cancer sample. 
As robust estimators we use the median and median absolute 
deviation (MAD; adjusted by a constant) 

OUTPUT: Three values seperated by a tab:
	1) median
	2) median absolute deviation
	3) number of paired-end-reads on which the estimates are based
"""

def convertPhredScore (phred_score):
	"""Converts phred score to probability"""
	return 10.0 ** (-float(phred_score) / float(10.0))

def abs_deviation(x,y):
	return abs(x - y)

def medianAbsoluteDeviation(array, median):
	"""Returns the MAD given the median."""
	compute_deviation = np.vectorize(abs_deviation)
	array = compute_deviation(array, median)
	return np.median(array)

def readIntegers(line):
	if line.strip() == '':
		return [] 
	l = []
	for elem in line.split():
		l.append(int(elem))
	return l 

def readFloats(line):
	if line.strip() == '':
		return [] 
	l = []
	for elem in line.split():
		l.append(float(elem))
	return l 

def exceedsThreshold(value, threshold):
	"""Returns true when value exceeds threshold. In case threshold = None, always returns True"""
	if threshold != None:
		if value < threshold: 
			return False
	return True

def isAutosome (chromosome_label):
	if chromosome_label[:3] == 'chr':
		chromosome_label = chromosome_label[3:]
	if chromosome_label in map(str, range(1,23)):
		return True
	return False
	

def main():
	parser = OptionParser(usage=usage)
	parser.add_option("--autosomes", action="store_true", dest="autosomes_only", default=False,
                      		help="Only insert size observations associated with autosomes are used.")
	parser.add_option("-q", action="store", type="float", dest="mapping_quality_threshold", default=None,
                      		help="Mapping quality threshold for the paired-end reads to be considered. (Default = None; i.e., no threshold).")
	parser.add_option("-v", action="store_true", dest="verbose", default=False,
                      		help="Verbose. Prints regularly how many variants have been processed.")
	(options, args) = parser.parse_args()

	if (len(args)!=1):
		parser.print_help()
		return 1

	observations_file = os.path.abspath(args[0])

	insert_sizes_healthy = [] 
	insert_sizes_tumour = [] 

	# walk through the variants: 
	n = 0
	with open(observations_file, 'r') as obs_file:
   		while True:
        		variant_data = list(islice(obs_file, 9)) # get 9 lines associated with one variant
        		if not variant_data:
           			break
			
			valid = True
			
			if options.autosomes_only: 
				# process the first line: 
				values 		= variant_data[0].split() 
				chromosome 	= values[1].strip()
				valid = isAutosome(chromosome) 

			
			if valid: 
				n += 1
				if options.verbose and n % 10000 == 0: 
					print('Having processed %d variants'%n, file=sys.stderr)

				# process the data 
				# HEALTHY SAMPLE 
				is_h 		= readIntegers(variant_data[1])	# insert sizes
				is_p_h	 	= readFloats(variant_data[2]) 	# insert sizes probabilities
		
				for insert_size, align_prob in zip(is_h, is_p_h):
					if exceedsThreshold(align_prob, options.mapping_quality_threshold):
						insert_sizes_healthy.append(insert_size) 	
			
				# TUMOUR SAMPLE 
				is_t  		= readIntegers(variant_data[5])	# insert sizes
				is_p_t  	= readFloats(variant_data[6])	# insert sizes probabilities

				for insert_size, align_prob in zip(is_t, is_p_t):
					if exceedsThreshold(align_prob, options.mapping_quality_threshold):
						insert_sizes_tumour.append(insert_size) 	


	insert_sizes_healthy 	= np.array(insert_sizes_healthy)
	insert_sizes_tumour 	= np.array(insert_sizes_tumour)

	median_healthy		= np.median(insert_sizes_healthy)	
	mad_healthy		= medianAbsoluteDeviation(insert_sizes_healthy, median_healthy)

	print("\n*** RESULTS ***\n")

	print("healthy sample")
	print("--------------")
	print("median: %lf\nMAD: %lf\n# observations: %d\n"%(median_healthy, mad_healthy * 1.4826, insert_sizes_healthy.size))

	median_tumour		= np.median(insert_sizes_tumour)	
	mad_tumour		= medianAbsoluteDeviation(insert_sizes_tumour, median_tumour)
	
	print("tumour sample")
	print("-------------")
	print("median: %lf\nMAD: %lf\n# observations: %d\n"%(median_tumour, mad_tumour * 1.4826, insert_sizes_tumour.size))

if __name__ == '__main__':
	sys.exit(main())
