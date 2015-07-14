#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
import os
import sys

from ResultsInstance import * 

__author__ = "Louis Dijkstra"

usage = """%prog <results-filename>

Summarizes results for deletions.

	<results-filename> 	File with the somatic mutation results. 
"""	

def main():

	parser = OptionParser(usage=usage)
	parser.add_option("--applyToAllChromosomes_WARNING_ONLY_IN_THE_END", action="store_true", dest="all_chromosomes", default=False,
						help="Summarize the results for all chromosomes (even and uneven), while normally only the even ones are used. ONLY USE IN THE END! This construction is used to avoid overfitting.")
	(options, args) = parser.parse_args()
	
	if (len(args)!=1):
		parser.print_help()
		return 1

	results, true_c_vafs = obtainResults(open(os.path.abspath(args[0]), 'r'), all_chromosomes = options.all_chromosomes, confidence_level = 0.0)

	for result in results:
		xh, yh, xc, yc = result.weightData()
		nxh, nyh, nxc, nyc = result.numberOfDataPoints()
		a, b = result.returnNumberSplitReads()
		call = result.call(p_value_threshold = 0.01, confidence_level = 0.0, confidence_interval = 0.95)
		width = result.ci_95_0[1] - result.ci_95_0[0]
		fi = result.determineObservedFisherInformation()
		print(result.obtainTruth(), '\t', result.length, '\t', xh, '\t', yh, '\t', xc, '\t', yc, '\t', nxh, '\t', nyh, '\t', nxc, '\t', nyc, '\t', a, '\t', b, '\t', width, '\t', call, '\t', result.ci_95_0[0], '\t', result.ci_95_0[1], '\t', fi, '\t', result.p_value)
					
if __name__ == '__main__':
	sys.exit(main())
