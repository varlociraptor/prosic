#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
from fractions import Fraction
import os
import sys
import math
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__))[:-3] + 'python')
from Indel import *

from matplotlib2tikz import save as tikz_save

__author__ = "Louis Dijkstra"

usage = """%prog <annotated-vcf-file>

	<annotated-vcf-file> 	A VCF file processes the somatic mutation caller.

Creates a box plot of the VAF estimates on the basis of the cancer BAM files. Considers both healthy and cancer cells. Used for assessing the VAF estimation quality. 
"""

def main():

	parser = OptionParser(usage=usage)
	parser.add_option("--even", action="store_true", dest="even_chromosomes", default=False,
						help="Summarize the results for all even chromosomes, while normally only the uneven ones are used. ONLY USE IN THE END! This construction is used to avoid overfitting.")
	parser.add_option("--final", action="store_true", dest="all_chromosomes", default=False,
						help="Summarize the results for all chromosomes, while normally only the uneven ones are used. ONLY USE IN THE END! This construction is used to avoid overfitting.")
	parser.add_option("--high", action="store_true", dest="high_precision", default=False,
						help="Uses high precision when matching detected deletions with the deletions in the true VCF file, i.e., the centerpoints must be less or exactly 20 bp away and the length must not differ more than 10 bp. This option overwrites the options '-d' and '-l' (!).")
	parser.add_option("-c", action="store", dest="confidence_level", default=0.0, type=float,
						help="Confidence level. (Default = 0.0; i.e., all)")
	parser.add_option("-d", action="store", dest="distance_threshold", default=100, type=int,
				  		help="Distance between centerpoints of deletions threshold. The distance must be smaller than this in order for the deletions to be deemed similar (in addition, the lengths must be similar, see option '-l'). (Default = 100 bp)")
	parser.add_option("-l", action="store", dest="length_threshold", default=100, type=int,
				  		help="Lengths of deletions threshold. The length difference must be smaller than this in order for the deletions to be deemed similar (in addition, the placement of the centerpoints must be similar, see option '-d'). (Default = 100 bp)")	
	parser.add_option("-o", action="store", dest="tikz_output", default=None, 
						help="Output is stored as a tikz file. (Default = not used)")
	parser.add_option("-t", action="store", dest="truth_summary_file", default=None, 
				  		help="The location of file that contains the ground truth. (Default = data/truth/truth.summary")
	parser.add_option("-x", action="store", dest="autosome", default=None, type=int, 
				  		help="Processes just one autosome.")
	(options, args) = parser.parse_args()
	
	if (len(args)!=1):
		parser.print_help()
		return 1

	# get truth file
	if options.truth_summary_file == None:
		options.truth_summary_file = os.path.abspath(os.path.dirname(__file__))[:-3] + 'data/truth/truth.summary'

	if options.high_precision:
		options.length_threshold 	= 10
		options.distance_threshold 	= 20

	list_of_autosomes = range(1,23,2) # only the uneven autosomes
	if options.even_chromosomes: 
		list_of_autosomes = range(2,23,2)
	if options.all_chromosomes:
		list_of_autosomes = range(1,23)
	if options.autosome != None:
		list_of_autosomes = [options.autosome]

	print('INPUT:\tfile with ground truth @', os.path.abspath(options.truth_summary_file), '\n\tcalls in VCF file @', os.path.abspath(args[0]), '\n')

	truth_file = open(options.truth_summary_file, 'r') # true indels 
	call_vcf_reader = vcf.Reader(open(os.path.abspath(args[0]))) # called indels

	




	results, true_c_vafs = obtainResults(open(os.path.abspath(args[0]), 'r'), all_chromosomes = options.all_chromosomes, confidence_level = options.confidence_level)

	# Collect the MLEs of the cancer VAFs
	c_vaf_mles = []
	for i in range(len(true_c_vafs)):
		c_vaf_mles.append([])
	
	for result in results:
		c_vaf_mles[true_c_vafs.index(result.true_c_vaf)].append(result.mle_c_vaf)

	plt.boxplot(c_vaf_mles, positions=true_c_vafs, widths = 0.02)
	xticks = [] 
	for x in true_c_vafs:
		 xticks.append(Fraction.from_float(x).limit_denominator().__str__())

	plt.xticks(true_c_vafs, xticks)

	#if options.coverage == 40:
	#	plt.xticks(range(1, 8), ['0', '1/12', '1/8', '1/3', '1/2', '2/3', '1'])
	#else:
	#	plt.xticks(range(1, len(true_c_vafs) + 1), true_c_vafs) # TODO make nicer
	plt.xlabel('true cancer VAF')
	plt.ylabel('cancer VAF MAP estimate')

	plt.plot(true_c_vafs, true_c_vafs, '*')

	x_min,x_max,y_min,y_max = plt.axis()
	plt.axis((-0.1, 1.1, y_min, y_max))

	#plt.show()	
	tikz_save( 'vaf_estimate_with_alpha.tikz', figureheight='4cm', figurewidth='6cm' )	
	
if __name__ == '__main__':
	sys.exit(main())
