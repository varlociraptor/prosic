#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
import os
import sys
import matplotlib.pyplot as plt
from collections import defaultdict
from itertools import islice

__author__ = "Louis Dijkstra"

usage = """%prog <observations-file> 

	<observations-file> 	observations file generated with
				'extract-observations.py' 

Plots the histogram of the insert distribution of both the healthy and cancer sample. 
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

	def print(self, output_stream):
		outputfile = open(output_stream, 'w') 
		if self.n == 0:
			return 0
		for i in range(self.min_value, self.max_value + 1):
			print("%d\t%d"%(i, self.count[i]), file=outputfile) 

	def returnData(self):
		x_values = range(self.min_value, self.max_value + 1)
		counts = []
		for x in x_values:
			counts.append(self.count[x])
		return x_values, counts 

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
	parser.add_option("--normalize", action="store_true", dest="normalization", default=False,  
				help="Normalizes the data.")
	parser.add_option("-k", action="store", dest="min_x", default=None, type=int,  
				help="Minimum x-value (Default: no minimum)")
	parser.add_option("-l", action="store", dest="max_x", default=None, type=int,   
				help="Maximum x-value (Default: no maximum)")
	parser.add_option("-o", action="store", dest="outputfilename", default=None,  
				help="Figure is stored under this name")	
	parser.add_option("-q", action="store", type="float", dest="mapping_quality_threshold", default=None,
                      		help="Mapping quality threshold for the paired-end reads to be considered. (Default = None; i.e., no threshold).")
	parser.add_option("-v", action="store_true", dest="verbose", default=False,
                      		help="Verbose. Prints regularly how many variants have been processed.")
	parser.add_option("-w", action="store", dest="width", default=1.0, type=float, 
				help="Bar width. (Default=1.0)")
	(options, args) = parser.parse_args()

	if (len(args)!=1):
		parser.print_help()
		return 1

	observations_file = os.path.abspath(args[0])

	health_is = CountData()
	tumour_is = CountData() 

	# walk through the observations file and collect the insert size observations 
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
						health_is.add(insert_size, 1) 	
			
				# TUMOUR SAMPLE 
				is_t  		= readIntegers(variant_data[5])	# insert sizes
				is_p_t  	= readFloats(variant_data[6])	# insert sizes probabilities

				for insert_size, align_prob in zip(is_t, is_p_t):
					if exceedsThreshold(align_prob, options.mapping_quality_threshold):
						tumour_is.add(insert_size, 1) 

	x_values_h, counts_h = health_is.returnData() 
	x_values_t, counts_t = tumour_is.returnData() 

	plt.xlabel("insert size")
	plt.ylabel("number of observations")
	
	if options.normalization: # normalize the count data 
		plt.ylabel("probability")
		total_obs_h = float(sum(counts_h))
		total_obs_t = float(sum(counts_t))
		for i, count in enumerate(counts_h):
			counts_h[i] = float(count) / total_obs_h 
		for i, count in enumerate(counts_t):
			counts_t[i] = float(count) / total_obs_t
		

	# Set the x-axis
	min_x = x_values_h[0] if (x_values_h[0] <= x_values_t[0]) else x_values_t[0]
	max_x = x_values_h[-1] if (x_values_h[-1] >= x_values_t[-1]) else x_values_t[-1]
	if options.min_x != None: 
		min_x = options.min_x 
	if options.max_x != None: 
		max_x = options.max_x 
	plt.xlim([min_x, max_x])

	data_h = plt.plot(x_values_h, counts_h, '.', color='blue', alpha=0.5,  label='healthy/control sample')
	data_t = plt.plot(x_values_t, counts_t, '.', color='red', alpha=0.5, label='tumour/case sample')
	plt.legend()
	#plt.bar(x_values, counts, align='center', width=options.width)

	if options.outputfilename != None: 
		plt.savefig(options.outputfilename, bbox_inches='tight')
	else:
		plt.show()
	

if __name__ == '__main__':
	sys.exit(main())
