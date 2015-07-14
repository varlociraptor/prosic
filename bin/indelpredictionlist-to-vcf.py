#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division
from optparse import OptionParser, OptionGroup
import sys
import os
from collections import defaultdict
from bisect import bisect_right

__author__ = "Tobias Marschall & Louis Dijkstra"

usage = """%prog [options] <indelprediction.list>

Outputs a VCF file in standard format (expected by PyVCF and 
DreamChallenge) given an indelprediction list outputed by 
Laser.
"""

def main():
	parser = OptionParser(usage=usage)
	parser.add_option("-m", action="store", dest="min_support", default=2, type=int,
					help='Minimal support (indels with less support are not printed). Default = 2')
	(options, args) = parser.parse_args()

	if (len(args) != 1):
		parser.print_help()
		sys.exit(1)

	print('##fileformat=VCFv4.1')
	print('##source=indelpredictionlist-to-vcf.py')
	print('##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">')
	print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO')
	for line in (s.strip() for s in open(args[0])):
		fields = line.split()
		chromosome, coord1, coord2, svtype, support, sequence = fields[0], int(fields[1]), int(fields[2]), fields[3].strip(), None, None
		if svtype == 'DEL':
			support = float(fields[4])
		elif svtype == 'INS':
			if len(fields) != 6:
				continue 
			sequence = fields[4]
			support = float(fields[5])
		else:
			assert False

		if support < options.min_support: # not enough support
			continue 

		if svtype == 'DEL':
			svlen = -(coord2-coord1+1)
			alt = '<DEL>'
			ref = '.'
			print('%s\t%d\t.\t%s\t%s\t.\tPASS\tSVTYPE=DEL;SVLEN=%d'%(chromosome, coord1, ref, alt, svlen))		
		elif svtype == 'INS':
			svlen = coord2
			alt = 'N' + sequence
			ref = '.' 
			print('%s\t%d\t.\t%s\t%s\t.\tPASS\tSVTYPE=INS;SVLEN=%d'%(chromosome, coord1, ref, alt, svlen))	
		else:
			assert False

if __name__ == '__main__':
	sys.exit(main())
