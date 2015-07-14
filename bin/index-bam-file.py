#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
import os
import sys
import pysam

__author__ = "Louis Dijkstra"

usage = """%prog <sorted-bam-file> 

	<sorted-bam-file> 	Sorted BAM file

Indexes a sorted BAM file sing PySAM (which employs SAMTools). 
One can create a sorted BAM file with the script 
'sort-bam-file.py'
"""

def main():

	parser = OptionParser(usage=usage)
	(options, args) = parser.parse_args()
	
	if (len(args)!=1):
		parser.print_help()
		return 1
		
	filename = args[0]
	
	if filename[-4:] != '.bam': 
		print("ERROR: file extension should be '.bam'. File is not indexed.")
		return 1 
		
	print("Indexing BAM file %s."%filename)
	pysam.index(filename)

if __name__ == '__main__':
	sys.exit(main())

