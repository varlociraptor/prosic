#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
from fractions import Fraction
import os
import sys
import math

import matplotlib.pyplot as plt
from ResultsInstance import *

__author__ = "Louis Dijkstra"

usage = """%prog <results-filename>

Creates a box plot of the VAF estimates on the basis of the cancer BAM files. Considers both healthy and cancer cells. Used for assessing the VAF estimation quality. 

	<results-filename> 	File with the somatic mutation results. 
"""

def main():

	parser = OptionParser(usage=usage)
	parser.add_option("--applyToAllChromosomes_WARNING_ONLY_IN_THE_END", action="store_true", dest="all_chromosomes", default=False,
						help="Summarize the results for all chromosomes (even and uneven), while normally only the even ones are used. ONLY USE IN THE END! This construction is used to avoid overfitting.")
	#parser.add_option("--bonf", action="store_true", dest="bonferroni", default=False,
	#					help="A Bonferroni correction is applied to the p-value threshold.")
	parser.add_option("-k", action="store", dest="length_start", default=None, type=int,
						help="Minimal length of the indels taken into account. (Default = unbounded)")
	parser.add_option("-l", action="store", dest="length_end", default=None, type=int,
						help="Maximal length of the indels taken into account. (Default = unbounded)")
	parser.add_option("-c", action="store", dest="confidence_level", default=0.95, type=float,
						help="Confidence level. (Default = 0.95)")
	#parser.add_option("-p", action="store", dest="p_value_threshold", default=0.05, type=float,
	#					help="P-value threshold used to make the call (somatic/not present). (Default = 0.05)")
	parser.add_option("-w", action="store", dest="max_width_ci", default=1.0, type=float,
						help="Threshold on the width of the 95% confidence interval. (Default = 1.0; all results are taken into account)")
	parser.add_option("-x", action="store", dest="coverage", default=40, type=int,
						help="Coverage of the cancer BAM file. (Default = 40)")
	
	(options, args) = parser.parse_args()
	
	if (len(args)!=1):
		parser.print_help()
		return 1

	results, true_c_vafs = obtainResults(open(os.path.abspath(args[0]), 'r'), all_chromosomes = options.all_chromosomes, confidence_level = options.confidence_level)

	# Collect the MLEs of the cancer VAFs
	c_vaf_mles = []
	for i in range(len(true_c_vafs)):
		c_vaf_mles.append([])
	
	#p_value_threshold = options.p_value_threshold
	#if options.bonferroni:
	#	p_value_threshold /= float(len(results))

	for result in results:
		mle, ci = result.returnCancerVAF_MLE()
		if ci[1] - ci[0] > options.max_width_ci: # if the confidence interval is too wide, ignore this particular result
			continue
		#if result.call(p_value_threshold = p_value_threshold, confidence_level = options.confidence_level, confidence_interval = 0.95) != 'somatic':
		#	continue
		c_vaf_mles[true_c_vafs.index(result.true_c_vaf)].append(mle)

	plt.boxplot(c_vaf_mles, positions=true_c_vafs, widths = 0.05)
	xticks = [] 
	for x in true_c_vafs:
		 xticks.append(Fraction.from_float(x).limit_denominator().__str__())


	plt.xticks(true_c_vafs, xticks)

	#if options.coverage == 40:
	#	plt.xticks(range(1, 8), ['0', '1/12', '1/8', '1/3', '1/2', '2/3', '1'])
	#else:
	#	plt.xticks(range(1, len(true_c_vafs) + 1), true_c_vafs) # TODO make nicer
	plt.xlabel('true cancer VAF')
	plt.ylabel('cancer VAF MLE')

	plt.plot(true_c_vafs, true_c_vafs, '*')

	x_min,x_max,y_min,y_max = plt.axis()
	plt.axis((-0.1, 1.1, y_min, y_max))

	
	plt.show()	

		
	
if __name__ == '__main__':
	sys.exit(main())
