#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
import os
import sys
import operator
from itertools import islice

__author__ = "Louis Dijkstra"

usage = """%prog [options] <observations-file> 

	<observations-file> 	original observations file 

(See 'extract-observations.py' for generating an observations file.)
Filters an observation file. Output is printed to standard output. 
See the options for the filter possibilities.
"""

def liesInInterval(value, min_value, max_value): 
	"""Returns True when the values lies within [min_value, max_value]. 
	   Otherwise False. In case, min_value or max_value is None, the interval
	   is considered open."""
	if min_value != None:
		if value < min_value: return False
	if max_value != None:
		if value > max_value: return False
	return True

def exceedsThreshold (value, threshold): 
	"""Returns whether value exceeds threshold. When threshold is None, 
	   it returns False"""
	if threshold != None: 
		return value > threshold 
	return False

def main():

	parser = OptionParser(usage=usage)
	parser.add_option("--deletions", action="store_true", dest="deletions", default=False,  
						help="Outputs only deletions.")
	parser.add_option("--insertions", action="store_true", dest="insertions", default=False,  
						help="Outputs only insertions")
	parser.add_option("-i", action="store", dest="del_is_threshold", default=None, type=int,
						help="Discards any insert size observations for any deletion that exceeds the given length. (Default = not applied)")
	parser.add_option("-I", action="store", dest="ins_is_threshold", default=None, type=int,
						help="Discards any insert size observations for any insertion that exceeds the given length. (Default = not applied)")
	parser.add_option("-k", action="store", dest="min_length", default=None, type=int,
						help="Outputs only indels larger or equal to this length. (Default = no restriction)")
	parser.add_option("-l", action="store", dest="max_length", default=None, type=int,
						help="Outputs only indels smaller or equal to this length. (Default = no restriction)")
	parser.add_option("-s", action="store", dest="del_split_threshold", default=None, type=int,
						help="Discards any split read observations for any deletion that exceeds the given length. (Default = not applied)")
	parser.add_option("-S", action="store", dest="ins_split_threshold", default=None, type=int,
						help="Discards any split read observations for any insertion that exceeds the given length. (Default = not applied)")
	parser.add_option("-x", action="store", dest="chromosome", default=None, 
						help="Outputs only variants on this chromosome.")
	(options, args) = parser.parse_args()
	
	if (len(args)!=1):
		parser.print_help()
		return 1

	observations_file = os.path.abspath(args[0])

	# walk through the variants: 
	with open(observations_file, 'r') as obs_file:
   		while True:
        		variant_data = list(islice(obs_file, 9)) # get 9 lines associated with one variant
        		if not variant_data:
           			break
			# process the first line: 
			values 		= variant_data[0].split() 
			variant_type 	= values[0].strip()
			chromosome 	= values[1].strip()
			position 	= int(values[2])
			length		= int(values[3])
			
			# check whether variant should be outputed or not	
			valid = liesInInterval(length, options.min_length, options.max_length) 	# length	
			if options.deletions:
				if variant_type != '-': valid = False	
			if options.insertions:
				if variant_type != '+': valid = False	
			if options.chromosome != None: 
				valid = valid and (chromosome == options.chromosome)

			if valid: # variant is valid and should be outputed:
				print(variant_data[0], end = '') # print data on variant 

				# check whether insert size observations or split observations should be neglected
				ignore_insert_size_obs, ignore_split_obs = False, False
				if variant_type == '-': # deletion
					ignore_insert_size_obs 	= exceedsThreshold(length, options.del_is_threshold)
					ignore_split_obs 	= exceedsThreshold(length, options.del_split_threshold)
				if variant_type == '+': # insertion
					ignore_insert_size_obs 	= exceedsThreshold(length, options.ins_is_threshold)
					ignore_split_obs 	= exceedsThreshold(length, options.ins_split_threshold)

				if ignore_insert_size_obs: print('\n') # 2 empty lines
				else: 
					print(variant_data[1], end = '')
					print(variant_data[2], end = '')
				
				if ignore_split_obs: print('\n') # 2 empty lines
				else: 
					print(variant_data[3], end = '')
					print(variant_data[4], end = '')

				if ignore_insert_size_obs: print('\n') # 2 empty lines
				else: 
					print(variant_data[5], end = '')
					print(variant_data[6], end = '')
				
				if ignore_split_obs: print('\n') # 2 empty lines
				else: 
					print(variant_data[7], end = '')
					print(variant_data[8], end = '')

if __name__ == '__main__':
	sys.exit(main())

