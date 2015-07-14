#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
import os
import sys
import pysam
import numpy as np
from collections import defaultdict

__author__ = "Louis Dijkstra"

usage = """%prog <bam-file> 

Returns robust estimates of the mean and standard deviation of the insert size
distribution. As robust estimators we use the median and median absolute 
deviation (MAD; adjusted by a constant) 

	<bam-file>	BAM file

OUTPUT: Three values seperated by a tab:
	1) median
	2) median absolute deviation
	3) number of paired-end-reads on which the estimates are based

NOTE: Assumes a constant read length for all reads in the BAM.
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

def determineReadLength(bam_file):
	"""Returns the read length of the reads in the BAM file. Assumes that all reads are of the same length."""
	for alignment in bam_file:
		if alignment.rlen is not None:
			return alignment.rlen
	return None

def returnAlignmentQuality(alignment1, alignment2):
	"""Returns the alignment quality (assumes independence)."""
	return (1.0 - convertPhredScore(alignment1.mapq)) * (1.0 - convertPhredScore(alignment2.mapq))

def main():
	parser = OptionParser(usage=usage)
	parser.add_option("-q", action="store", type="float", dest="mapping_quality_threshold", default=None,
                      		help="Mapping quality threshold for the paired-end reads to be considered. (Default = None; i.e., no threshold).")
	parser.add_option("-v", action="store_true", dest="verbose", default=False,
                      		help="Verbose. Prints regularly how many alignments have been processed.")
	(options, args) = parser.parse_args()

	if (len(args)!=1):
		parser.print_help()
		return 1

	input_file = pysam.Samfile(args[0], "rb")
	read_lengths = 2 * determineReadLength(input_file) # total length of paired-end read	

	insert_sizes = [] # insert sizes of paired-end read are stored here

	n = 0 
	if options.mapping_quality_threshold == None: # all alignments are taken into account
		for alignment in input_file:
			n += 1
			if options.verbose and n % 100000 == 0: 
				print('Having processed %d alignments'%n, file=sys.stderr)
			if alignment.isize - read_lengths <= 0:
				continue
			insert_sizes.append(alignment.isize - read_lengths)
	else: # paired alignments for which the alignment probability is high enough are taken into account
		alignment_dict = defaultdict(list)
		# pair the alignments, when possible
		for alignment in input_file:
			n += 1
			if options.verbose and n % 100000 == 0: 
				print('Having processed %d alignments'%n, file=sys.stderr)
			if alignment.isize == 0: # alignment is not mapped
				continue
			alignment_dict[alignment.qname].append(alignment)
		
		# walk through all paired alignments
		for qname, alignments in alignment_dict.iteritems():
			if len(alignments) == 2: # paired-end read
				if returnAlignmentQuality(alignments[0], alignments[1]) >= options.mapping_quality_threshold:
					insert_sizes.append(abs(alignments[0].isize) - read_lengths)

	insert_sizes 	= np.array(insert_sizes)
	median 		= np.median(insert_sizes)	
	mad		= medianAbsoluteDeviation(insert_sizes, median)

	print(int(median), '\t', mad * 1.4826, '\t', insert_sizes.size)

if __name__ == '__main__':
	sys.exit(main())
