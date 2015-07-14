#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
import os
import sys
import math

from summarize_results import *
from ResultsInstance import *
from ternary_classification import * 

__author__ = "Louis Dijkstra"

usage = """%prog <results-filename>

Plots the calling results (SOMATIC/GERMLINE/NOT PRESENT). The length ranges can be changed in the code (LENGTH_RANGES)

	<results-filename> 	File with the somatic mutation results. 
"""

LENGTH_RANGES = [[10, None], [10,29], [30,49], [50, 69], [70, 100]]

def main():
	parser = OptionParser(usage=usage)
	parser.add_option("--applyToAllChromosomes_WARNING_ONLY_IN_THE_END", action="store_true", dest="all_chromosomes", default=False,
						help="The results for all autosomes are used, while normally only the first chromosome is used. ONLY USE IN THE END! This construction is used to avoid overfitting.")
	parser.add_option("--bonf", action="store_true", dest="bonferroni", default=False,
						help="A Bonferroni correction is applied to the p-value threshold.")
	parser.add_option("--use-ci", action="store_true", dest="use_ci", default=False,
						help="The 95% confidence interval is used to make the call, i.e., when 0 is in the 95% CI, we call it 'not present', otherwise 'somatic'.")
	parser.add_option("-c", action="store", dest="confidence_level", default=0.95, type=float,
						help="Confidence level. (Default = 0.95)")
	parser.add_option("-d", action="store", dest="width", default=0.9, type=float,
						help="Bar width. (Default = 0.9)")
	parser.add_option("-i", action="store", dest="internal_segment_threshold", default=0.0, type=float, 
						help="Number of internal segment alignments from the cancer sample. (Default = 0.0, i.e., no threshold).")
	parser.add_option("-p", action="store", dest="p_value_threshold", default=0.05, type=float, 
						help="P-value threshold used for deciding whether an indel is somatic/not present. (Default = 0.05).")
	parser.add_option("-s", action="store", dest="split_read_threshold", default=0.0, type=float, 
						help="Number of split reads from the cancer sample that support the presence of the deletion. (Default = 0.0, i.e., no threshold).")
	parser.add_option("-w", action="store", dest="max_width_ci", default=1.0, type=float,
						help="Threshold on the width of the confidence interval. (Default = 1.0; all results are taken into account)")
	(options, args) = parser.parse_args()
	
	if (len(args)!=1):
		parser.print_help()
		return 1

	results, true_c_vafs = obtainResults(open(os.path.abspath(args[0]), 'r'), all_chromosomes = options.all_chromosomes, confidence_level = 0.0, threshold_internal_segment_evidences = options.internal_segment_threshold, threshold_split_read_evidences = options.split_read_threshold)
	
	p_value_threshold = options.p_value_threshold 
	if options.bonferroni:
		p_value_threshold /= float(len(results))

	somatic, germline, not_present = [], [], []  
	for length_range in LENGTH_RANGES:

		percentage_correctly_called, n_test_cases_used, total, table = testCallingPerformance(results, length_range, options.confidence_level, [None, None], 0.95, p_value_threshold, options.max_width_ci, use_ci = options.use_ci)	
		table = normalizeTable(table)
		somatic.append(table[0])
		germline.append(table[1])
		not_present.append(table[2])
	
	plotTernaryClassification(somatic, germline, not_present, LENGTH_RANGES, class_names = ['somatic', 'germline', 'not present'], width = options.width)			
	
if __name__ == '__main__':
	sys.exit(main())
