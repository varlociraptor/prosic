#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
import os
import sys
import pysam



__author__ = "Louis Dijkstra"

usage = """%prog <bam-file> 

	Sorts and indexes a BAM file. New BAM file is stored with the extension
	'sorted.bam'

	<bam-file>	BAM file that needs to be sorted and indexed
"""


def main():


	parser = OptionParser(usage=usage)
	(options, args) = parser.parse_args()

	if (len(args)!=1):
		parser.print_help()
		return 1

	bam_filename = os.path.abspath(args[0])

	if bam_filename[-4:] != ".bam":
		print("ERROR: file should have the extension .bam\nFile not sorted & indexed.")
		return 1

	sorted_bamfilename = bam_filename[:-4] + ".sorted"
	
	print("Sorting BAM file ", bam_filename)
	pysam.sort(bam_filename, sorted_bamfilename)
	print("Indexing the now sorted BAM file:", sorted_bamfilename + ".bam")
	pysam.index(sorted_bamfilename + ".bam")
	print("Done")

if __name__ == '__main__':
	sys.exit(main())
