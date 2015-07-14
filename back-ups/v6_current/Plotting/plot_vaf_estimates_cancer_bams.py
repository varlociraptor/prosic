#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
from fractions import Fraction
import os
import sys
import math
from matplotlib2tikz import save as tikz_save

import matplotlib.pyplot as plt
from ResultsInstance import *

__author__ = "Louis Dijkstra"

usage = """%prog <results-filename>

Creates a box plot of the VAF estimates on the basis of the cancer BAM files. Considers both healthy and cancer cells. Used for assessing the VAF estimation quality. 

	<results-filename> 	File with the somatic mutation results. 
"""

def determineCoverage(filename):
	"""Returns the coverage on the basis of the filename."""
	if '40x' in filename:
		return 40
	if '80x' in filename:
		return 80
	print('ERROR: Coverage is unknown. Filename should contain the symbols 40x or 80x.')
	sys.exit()

def main():

	parser = OptionParser(usage=usage)
	parser.add_option("--applyToAllChromosomes_WARNING_ONLY_IN_THE_END", action="store_true", dest="all_chromosomes", default=False,
						help="Summarize the results for all chromosomes (even and uneven), while normally only the even ones are used. ONLY USE IN THE END! This construction is used to avoid overfitting.")
	parser.add_option("-k", action="store", dest="length_start", default=None, type=int,
						help="Minimal length of the indels taken into account. (Default = unbounded)")
	parser.add_option("-l", action="store", dest="length_end", default=None, type=int,
						help="Maximal length of the indels taken into account. (Default = unbounded)")
	parser.add_option("-c", action="store", dest="confidence_level", default=0.95, type=float,
						help="Confidence level. (Default = 0.95)")
	parser.add_option("-w", action="store", dest="max_width_ci", default=1.0, type=float,
						help="Threshold on the width of the 95% confidence interval. (Default = 1.0; all results are taken into account)")
	(options, args) = parser.parse_args()
	
	if (len(args)!=1):
		parser.print_help()
		return 1

	filename = os.path.abspath(args[0])
	coverage = determineCoverage(filename)

	results, true_c_vafs = obtainResults(open(filename, 'r'), all_chromosomes = options.all_chromosomes, confidence_level = options.confidence_level)

	true_vafs = [h/16.0 for h in range(16+1)] # possible VAFs one can find in the cancer BAM when the coverage is 40x
	alpha = .25 # level of impurity for the cancer 40x and 80x BAM
	if coverage == 80:
		true_vafs = [h/32.0 for h in range(32+1)]  # possible VAFs one can find in the cancer BAM when the coverage is 80x

	vaf_mles = [] 
	for i in range(len(true_vafs)):
		vaf_mles.append([])

	for result in results:
		true_vaf = round(result.true_h_vaf*alpha + (1 - alpha)*result.true_c_vaf, 5)
		vaf_mles[true_vafs.index(true_vaf)].append(result.mle_c_vaf)

	plt.boxplot(vaf_mles, positions=true_vafs, widths = 0.02)
	
	xticks = [] 
	for x in true_c_vafs:
		 xticks.append(Fraction.from_float(x).limit_denominator().__str__())

	plt.xticks(true_c_vafs, xticks)
	plt.xlabel('true VAF')
	plt.ylabel('MLE')

	dots = [] 
	for result in results:
		true_vaf = round(result.true_h_vaf*alpha + (1 - alpha)*result.true_c_vaf, 5)
		if true_vaf not in dots:
			dots.append(true_vaf)

	plt.plot(dots, dots, '*')
	x_min,x_max,y_min,y_max = plt.axis()
	plt.axis((-0.1, 1.1, y_min, y_max))

	#plt.show()			
	tikz_save( 'vaf_estimate.tikz', figureheight='4cm', figurewidth='6cm' )
	
if __name__ == '__main__':
	sys.exit(main())
