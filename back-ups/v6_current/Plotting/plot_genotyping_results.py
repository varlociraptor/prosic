#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
import os
import sys
import math

from ResultsInstance import *
from TernaryClassification import * 

__author__ = "Louis Dijkstra"

usage = """%prog <results-filename>

Plots the genotyping results for the Control BAM. The length ranges can be changed in the code (LENGTH_RANGES)

	<results-filename> 	File with the somatic mutation results. 
"""

LENGTH_RANGES = [[10, None], [10,29], [30,49], [50,69], [70,100]]


def returnIndex (chromosome):
	if len(chromosome) > 3 and chromosome[:3] == 'chr':
		chromosome = chromosome[3:]
	if chromosome == 'X' or chromosome == 'x':
		return 22
	if chromosome == 'Y' or chromosome == 'y':
		return 23	
	return int(chromosome)

def main():

	parser = OptionParser(usage=usage)
	parser.add_option("--applyToAllChromosomes_WARNING_ONLY_IN_THE_END", action="store_true", dest="all_chromosomes", default=False,
						help="The results for all autosomes are used, while normally only the first chromosome is used. ONLY USE IN THE END! This construction is used to avoid overfitting.")
	parser.add_option("-c", action="store", dest="confidence_level", default=0.0, type=float,
						help="Confidence level. (Default = 0.0)")
	parser.add_option("-w", action="store", dest="width", default=0.9, type=float,
						help="Bar width. (Default = 0.9)")
	(options, args) = parser.parse_args()
	
	if (len(args)!=1):
		parser.print_help()
		return 1

	results, true_c_vafs = obtainResults(open(os.path.abspath(args[0]), 'r'), all_chromosomes = options.all_chromosomes, confidence_level = options.confidence_level)

	# TODO remove later
	list_of_autosomes = range(1,2) # TODO change back!
	results_new = []
	for result in results:
		if returnIndex(result.chromosome) in list_of_autosomes:
			results_new.append(result)
	results = results_new

	not_present, heterozygous, homozygous = [], [], []  
	for length_range in LENGTH_RANGES:
		total = 0.0
		table = [[0,0,0],[0,0,0],[0,0,0]]
		for result in results: 
			if not result.fallsInLengthRange(lower=length_range[0], upper=length_range[1]):
				continue 
			genotype = result.genotypeControl(confidence_level = options.confidence_level) 
			if genotype is None:
				continue
			row, column = 0,0 
			if genotype == 0.5: 		row = 1
			if genotype == 1.0: 		row = 2
			if result.true_h_vaf == 0.5: 	column = 1
			if result.true_h_vaf == 1.0:	column = 2
			table[row][column] += 1
		table = normalizeTable(table)
		not_present.append(table[0])
		heterozygous.append(table[1])
		homozygous.append(table[2])
	
	plotTernaryClassification (not_present, heterozygous, homozygous, LENGTH_RANGES, class_names = ['not present', 'heterozygous', 'homozygous'], width = options.width)	
	
if __name__ == '__main__':
	sys.exit(main())
