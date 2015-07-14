#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division
from optparse import OptionParser, OptionGroup
import sys
import os
from collections import defaultdict
from bisect import bisect_right

__author__ = "Tobias Marschall"

usage = """%prog [options] <deletions.list>"""

def main():
	parser = OptionParser(usage=usage)
	parser.add_option("-m", action="store", dest="min_support", default=2, type=int,
					help='Minimal support (indels with less support are not printed).')
	(options, args) = parser.parse_args()
	if (len(args) != 1):
		parser.print_help()
		sys.exit(1)
	print('##fileformat=VCFv4.1')
	print('##source=deletionlist-to-vcf.py')
	print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tdefault')
	for line in (s.strip() for s in open(args[0])):
		fields = line.split()
	#	assert len(fields) == 5
		chromosome, coord1, coord2, svtype, support = fields[0], int(fields[1]), int(fields[2]), fields[3], fields[4]
#		chromosome, coord1, coord2, svtype = fields[0], int(fields[1]), int(fields[2]), fields[3]
		if float(support) < options.min_support:
			continue
                support = 1
		if svtype == 'DEL':
			print(chromosome, coord1, '.', '.', '<DEL>', '.', 'PASS', 'SVTYPE=DEL;SVLEN=%d'%(-(coord2-coord1+1)), 'GT:AD', '1/.:%d'%support)
		elif svtype == 'INS':
			print(chromosome, coord1, '.', '.', '<DEL>', '.', 'PASS', 'SVTYPE=INS;SVLEN=%d'%coord2, 'GT:AD', '1/.:%d'%support)
		else:
			assert False

if __name__ == '__main__':
	sys.exit(main())
