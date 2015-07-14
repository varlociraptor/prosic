#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
import os
import sys
import pysam

__author__ = "Louis Dijkstra"

usage = """%prog <bam-file> <output-bam-file>

Sorts and indexed the BAM file using pysam. 
 
	<bam-file> 		original BAM file
	<output-bam-file>	name where the sorted and indexed BAM will be stored	
"""

def main():

	parser = OptionParser(usage=usage)
	(options, args) = parser.parse_args()
	
	if (len(args)!=2):
		parser.print_help()
		return 1

	pysam.sort(args[0], args[1][:-3])
	pysam.index(args[1])

if __name__ == '__main__':
	sys.exit(main())

