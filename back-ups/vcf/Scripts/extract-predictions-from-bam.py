#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function, division
from optparse import OptionParser
import sys
import os
from pysam import Samfile, AlignedRead
from collections import defaultdict

__author__ = "Tobias Marschall"

usage = """%prog [options]

Reads BAM from stdin and outputs insertion and deletion predictions present in the alignments.
Coordinates are 0-based and (in case of deletions) inclusive."""

def main():
	parser = OptionParser(usage=usage)
	parser.add_option("-s", action="store_true", dest="sam_input", default=False,
				  help="Input is in SAM format instead of BAM format")
	parser.add_option("-m", action="store", dest="min_length", default=10, type=int,
				  help="Minimal size of events to be extracted (default=10).")
	parser.add_option("-n", action="store", dest="readnames_file", default=None,
				  help="Name of file to which readnames of split-reads are to be written.")
	(options, args) = parser.parse_args()
	if (len(args) != 0) or (os.isatty(0)):
		parser.print_help()
		sys.exit(1)
	readnames_file = None
	if options.readnames_file != None:
		readnames_file = open(options.readnames_file, 'w')
	inputfile = Samfile('-', 'r' if options.sam_input else 'rb')
	# maps tuples (chr,start,end) to number of times observed
	deletions = defaultdict(int)
	# maps tuples (chr,start,length) to number of times observed
	insertions = defaultdict(int)
	n = 0
	for read_aln in inputfile:
	        n += 1
	        if n % 1000000 == 0:
	                print('Having processed %d alignments'%n, file=sys.stderr)
		if read_aln.is_unmapped: continue
		if read_aln.is_secondary: continue
		ref_pos = read_aln.pos
		is_split_read = False
		for event, length in read_aln.cigar:
			if event in [0,7,8]: # Match/mismatch
				ref_pos += length
			elif event == 1: # Insertion
				if length >= options.min_length:
					insertions[(read_aln.tid, ref_pos, length)] += 1
					is_split_read = True
					#print(inputfile.getrname(read_aln.tid), ref_pos+1, length, 'INS')
			elif event == 2: # Deletion
				if length >= options.min_length:
					deletions[(read_aln.tid, ref_pos, ref_pos+length)] += 1
					is_split_read = True
					#print(inputfile.getrname(read_aln.tid), ref_pos+1, ref_pos+length, 'DEL')
				ref_pos += length
			elif event == 3: # Refskip
				ref_pos += length
			elif event == 4: # Softclip
				pass
			elif event == 5: # Hardclip
				pass
			else: 
				assert False
		if is_split_read and (readnames_file != None): 
			print(read_aln.qname, file=readnames_file)
	l = [(chr_id, start, end, count) for (chr_id, start, end), count in deletions.iteritems()]
	l.sort()
	for chr_id, start, end, count in l:
		print(inputfile.getrname(chr_id), start, end-1, 'DEL', count)
	l = [(chr_id, start, length, count) for (chr_id, start, length), count in insertions.iteritems()]
	l.sort()
	for chr_id, start, length, count in l:
		print(inputfile.getrname(chr_id), start, length, 'INS', count)
	#l = [(chr_id, start, end, count) for count,(chr_id, start, end) in enumerate ]


if __name__ == '__main__':
	sys.exit(main())
