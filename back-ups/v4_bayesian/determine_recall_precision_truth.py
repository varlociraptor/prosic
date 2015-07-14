#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
import vcf
import os
import sys
import math

from Indel import *
from Alignments import * 
from ResultsInstance import *

__author__ = "Louis Dijkstra"

usage = """%prog options <results-filename>

Determines the recall and precision (as defined in the MATE-CLEVER paper) for the deletions in the 
true VCF file. Only deletions on the autosomes are taken into account. The length ranges can be changed manually
in the code. 

	<results-filename> 	File with the somatic mutation results.
"""

LENGTH_RANGES 		= [[None, None], [10,29], [30,49], [50,69], [70,99], [100, 249], [250, None]]

VCF_TRUTH_FILENAME = "/data1/structvar/vaf-experiments/vcf/cancerclones.vcf.gz" # VCF file that contains the ground truth
VCF_TRUTH_READER = vcf.Reader(open(VCF_TRUTH_FILENAME)) # used for reading the VCF file 

list_true_deletions = [[] for i in range(24)] # 24 empty lists for the 22 autosomes and the two sex-chromosomes X and Y   

def returnIntervalString (interval):
	if interval[0] is None and interval[1] is not None:
		return "at most " + str(interval[1])
	elif interval[0] is not None and interval[1] is None:
		return "at least " + str(interval[0])
	elif interval[0] is None and interval[1] is None:
		return "unbounded"
	else:
		return "[" + str(interval[0]) + ','  + str(interval[1]) + "]"

def returnIndex (chromosome):
	if len(chromosome) > 3 and chromosome[:3] == 'chr':
		chromosome = chromosome[3:]
	if chromosome == 'X' or chromosome == 'x':
		return 22
	if chromosome == 'Y' or chromosome == 'y':
		return 23	
	return int(chromosome)

def obtainListOfTrueDeletions (min_length): 
	"""Fills the global variable list_true_deletions with the true deletions from the true VCF file"""
	for true_vcf_record in VCF_TRUTH_READER:
		if not isDeletion(true_vcf_record): # not a deletion
			continue 
		if returnIndelLength(true_vcf_record) < min_length: # too small
			continue 
		list_true_deletions[returnIndex(true_vcf_record.CHROM)].append(true_vcf_record) 

def deletionsSimilar(deletion, result, difference_length_threshold = 100, difference_centerpoint_threshold = 100):
	"""Returns true when deletion and result are similar, otherwise false."""
	# determine centerpoints
	centerpoints_other_del = [] # centerpoints of the other deletion
	if result.length % 2 == 0: 
		centerpoints_other_del = [result.position + 1 + result.length / 2 - 1, result.position + 1 + result.length / 2]
	else: 
		centerpoints_other_del = [result.position + 1 + result.length / 2]
	if returnMinimumDifference(deletion.centerpoints, centerpoints_other_del) > difference_centerpoint_threshold:
		return False
	if abs(deletion.length - result.length) > difference_length_threshold:
		return False
	return True

def main():

	parser = OptionParser(usage=usage)
	parser.add_option("--low", action="store_true", dest="low_precision", default=False,
						help="Uses low precision when matching detected indels with the indels in the true VCF file, i.e., the centerpoints must be less or exactly 100 bp away and the length must not differ more than 100 bp. This option overwrites the options '-d' and '-l' (!).")
	parser.add_option("-d", action="store", dest="distance_threshold", default=100, type=int,
				  		help="Distance between centerpoints of deletions threshold. The distance must be smaller than this in order for the deletions to be deemed similar (in addition, the lengths must be similar, see option '-l'). (Default = 100 bp)")
	parser.add_option("-i", action="store", dest="internal_segment_threshold", default=0.0, type=float, 
						help="Number of internal segment alignments from the cancer sample. (Default = 0.0, i.e., no threshold).")
	parser.add_option("-l", action="store", dest="length_threshold", default=100, type=int,
				  		help="Lengths of deletions threshold. The length difference must be smaller than this in order for the deletions to be deemed similar (in addition, the placement of the centerpoints must be similar, see option '-d'). (Default = 100 bp)")
	parser.add_option("-m", action="store", dest="min_length", default=10, type=int,
				  		help="Minimal length of an indel to be considered. (Default = 10 bp)")
	parser.add_option("-s", action="store", dest="split_read_threshold", default=0.0, type=float, 
						help="Number of split reads from the cancer sample that support the presence of the deletion. (Default = 0.0, i.e., no threshold).")
	(options, args) = parser.parse_args()

	if (len(args)!=1):
		parser.print_help()
		return 1

	if options.low_precision:
		options.length_threshold = 100
		options.distance_threshold = 100

	print('Reading true VCF file @ /data1/structvar/vaf-experiments/vcf/cancerclones.vcf.gz')
	obtainListOfTrueDeletions(options.min_length) # read in all the true deletions from the true VCF file
	print('Done with reading in the true VCF.')
	print('Reading in the results @ ' + str(os.path.abspath(args[0])))
	results, true_c_vafs = obtainResults(open(os.path.abspath(args[0]), 'r'), all_chromosomes = True, confidence_level = 0.0, threshold_internal_segment_evidences = options.internal_segment_threshold, threshold_split_read_evidences = options.split_read_threshold)
	print('Done with reading in the results.')

	n_true_annotations, n_predictions = [0 for i in range(len(LENGTH_RANGES))], [0 for i in range(len(LENGTH_RANGES))]
	n_true_positives = [0 for i in range(len(LENGTH_RANGES))] # initialize 

	results_per_chromosome = [[] for i in range(24)]
	for result in results:
		results_per_chromosome[returnIndex(result.chromosome)].append(result)
		for i in range(len(LENGTH_RANGES)):
			length_range = LENGTH_RANGES[i]
			if result.fallsInLengthRange(length_range[0], length_range[1]):
				n_predictions[i] += 1

	# walk through all autosomes (sex chromosomes are not taken into account)
	for autosome in range(1,23): 
		print('Processing autosome', autosome)
		# walk through all true deletions and determine they have been called or not
		n_true_positives_chr = 0 
		for vcf_record in list_true_deletions[autosome]:
			true_deletion = Deletion(vcf_record)
			called_correctly = False
			for result in results_per_chromosome[autosome]:
				if deletionsSimilar(true_deletion, result, difference_length_threshold = options.length_threshold, difference_centerpoint_threshold = options.distance_threshold):
					called_correctly = True
		
			for i in range(len(LENGTH_RANGES)):
				length_range = LENGTH_RANGES[i]
				if true_deletion.fallsInLengthRange(length_range[0], length_range[1]):
					n_true_annotations[i] += 1
					if called_correctly:
						n_true_positives[i] += 1
	
	
	for i in range(len(LENGTH_RANGES)):
		print("LENGTH RANGE: ", returnIntervalString(LENGTH_RANGES[i]))	
		print('total number of true deletions in the VCF file: ', n_true_annotations[i])
		print('total number of predictions made: ', n_predictions[i])
		print('true positives: ', n_true_positives[i])
		print('recall: ', float(n_true_positives[i]) / float(n_true_annotations[i]))
		print('precision: ', float(n_true_positives[i]) / float(n_predictions[i]))		
		print('')
		
if __name__ == '__main__':
	sys.exit(main())




