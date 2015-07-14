#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function, division
from optparse import OptionParser
import sys
import os
from pysam import Samfile
import bisect

__author__ = "Tobias Marschall"

usage = """%prog [options]

Reads BAM format from stdin and outputs a list of plausible alignments
(one per line) of read pairs to stdout in the following format:
<read-name> <chromosome1> <start1> <end1> <strand1> <chromosome2> <start2> <end2> <strand2> <aln-pair-prob> <aln-pair-prob-norm>
Where <aln-pair-prob> is the (estimated) probability that this alignment pair is 
correct. The probabilities given in <aln-pair-prob-norm> are normalized such
that they sum up to one for each read pair.

<startX> and <endX> coordinates are 1-based and inclusive (i.e. closed intervals)
and give the region on the reference the respective read was aligned to.

The probability estimates are based on the alignment score (AS) tags in the 
BAM input. Therefore, these tags must be present.

Unmapped reads are ignored. Alignment pairs for which both alignments have the same 
orientation (strandedness) or map to different chromosomes are not reported. 

IMPORTANT: Assumes that, in the BAM input, alignments belonging to the same read
           are grouped, i.e. are given in subsequent lines.

NOTE: When using BWA, multiple alignment locations are by default stored in XA tags
      in a single line (i.e. one line per read rather than per alignment). Can be 
      converted using the script xa2multi.pl that comes with BWA."""

skipped_duplicates = 0

def length_on_reference(aln_read):
	"""Returns the number of letters in the references that are covered by the given read alignment."""
	# CIGAR codes:
	# BAM_CMATCH     = 0 (M)
	# BAM_CINS       = 1 (I)
	# BAM_CDEL       = 2 (D)
	# BAM_CREF_SKIP  = 3 (N)
	# BAM_CSOFT_CLIP = 4 (S)
	# BAM_CHARD_CLIP = 5 (H)
	# BAM_CPAD       = 6 (P)
	if aln_read.is_unmapped: return 0
	return aln_read.rlen + sum(l for id,l in aln_read.cigar if id == 2) - sum(l for id,l in aln_read.cigar if id == 1)

def compute_posteriors(alignments):
	"""Attaches probabilities corresponding with the 'AS' scores
	of the alignments."""
	scores = []
	for read_aln in alignments:
		asflag = False
		for tag in read_aln.tags:
			if tag[0] == 'AS':
				scores.append(tag[1])
				asflag = True
				break
		if not asflag:
			scores.append(None)
	probsum = 0.0
	probs = []
	for score in scores:
		if score == None:
			prior = 0.0
		else:
			prior = 10.0**(score/-10.0)
		probs.append(prior)
		probsum += prior
	if probsum > 0.0:
		probs = [x/probsum for x in probs]
	else:
		probs = [1.0/len(alignments) for x in probs]
	return probs

def process_read_allpairs(read_name, refnames, alignments1, alignments2, insert_length_distribution):
	"""Given an output file (a plain file object) and a list of alignments belonging
	to the same read pair, outputs one line per alignment pair in the format:
	"<read-name> <chromosome1> <pos1> <strand1> <chromosome2> <pos2> <strand2> <aln-pair-prob> <aln-pair-prob-norm>"."""
	posteriors1 = compute_posteriors(alignments1)
	posteriors2 = compute_posteriors(alignments2)
	assert len(posteriors1) == len(alignments1)
	assert len(posteriors2) == len(posteriors2)
	# contains a tuple (read_aln1, read_aln2, pair-prob) for every pair
	# that should be printed
	all_pairs = []
	psum = 0.0
	for read_aln1, posterior1 in zip(alignments1,posteriors1):
		assert not read_aln1.is_unmapped
		for read_aln2, posterior2 in zip(alignments2,posteriors2):
			assert not read_aln2.is_unmapped
			# skip if strandedness is the same
			if read_aln1.is_reverse == read_aln2.is_reverse: continue
			# skip if both reads map to different chromosomes
			if read_aln1.rname != read_aln2.rname: continue
			p = posterior1*posterior2
			if read_aln1.pos <= read_aln2.pos:
				ra1, ra2 = read_aln1, read_aln2
			else:
				ra1, ra2 = read_aln2, read_aln1
			if insert_length_distribution != None:
				insert_length = ra2.pos - ra1.pos + length_on_reference(ra1)
				p *= insert_length_distribution.probability(insert_length)
			psum += p
			all_pairs.append((ra1, ra2, p))
	if psum == 0.0:
		return
	for read_aln1, read_aln2, pair_posterior in all_pairs:
		print(read_name, 
			refnames[read_aln1.tid], read_aln1.pos+1, read_aln1.pos + length_on_reference(read_aln1), '-' if read_aln1.is_reverse else '+', 
			refnames[read_aln2.tid], read_aln2.pos+1, read_aln2.pos + length_on_reference(read_aln2), '-' if read_aln2.is_reverse else '+', 
			pair_posterior, pair_posterior/psum
		)

class InsertLengthDistribution:
	def __init__(self, filename):
		"""Reads a distribution from a file whose lines have the following format: <start> <end> <probability>,
		where the last column must add up to one."""
		values = [(int(start),int(end),float(probability)) for start,end,probability in (s.split() for s in file(filename))]
		# do some checks
		last = None
		psum = 0.0
		for start, end, probability in values:
			assert (last==None) or (start==last+1)
			assert start <= end
			last = end
			psum += probability
		assert abs(1.0 - psum) < 1e-10
		self.starts = [start for start,end,probability in values]
		self.probabilities = [probability/(end-start+1) for start,end,probability in values]
		self.min = values[0][0]
		self.max = values[-1][1]
	def __str__(self):
		return str(zip(self.starts,self.probabilities))
	def probability(self, insert_length):
		if insert_length < self.min: return 0.0
		if insert_length > self.max: return 0.0
		i = bisect.bisect_right(self.starts, insert_length) - 1
		return self.probabilities[i]

def main():
	print("This script is obsolete. Use the C++ version.", file=sys.cerr)
	sys.exit(1)

	global skipped_duplicates
	parser = OptionParser(usage=usage)
	parser.add_option("-s", action="store_true", dest="sam_input", default=False,
					  help="Input is in SAM format instead of BAM format")
	parser.add_option("-V", action="store_true", dest="verify_sortedness", default=False,
					  help="Verifies that alignments for the same read are indeed grouped together (at the expense of using more RAM)")
	parser.add_option("-X", action="store_true", dest="skip_non_xa", default=False,
					  help="Skip reads for which other alignments exist (i.e. X0+X1>1), but no XA tag is present.")
	parser.add_option("-i", action="store", dest="insert_length_histogram",
					  help="Take insert length distribution into account. The given file must contain a histogram as produced by \"insert-length-histogram.py -B\".")
	(options, args) = parser.parse_args()
	if (len(args)!=0) or (os.isatty(0)):
		parser.print_help()
		sys.exit(1)
	inputfile = Samfile('-', 'r' if options.sam_input else 'rb')
	if options.verify_sortedness:
		read_names = set()
	insert_length_distribution = None
	if options.insert_length_histogram != None:
		insert_length_distribution = InsertLengthDistribution(options.insert_length_histogram)
	last_read_name = None
	# alignments of the currently processed read pair
	alignments1, alignments2 = [], []
	alignments1_starts, alignments2_starts = set(), set()
	counter = 0
	skipped_by_xa = 0
	for read_aln in inputfile:
		counter += 1
		if counter % 1000000 == 0:
			print("having processed %d read alignments" % (counter), file=sys.stderr)
		assert read_aln.is_read1 != read_aln.is_read2
		if read_aln.is_unmapped: continue
		if options.skip_non_xa:
			n = 0
			has_xa = False
			for tag, value in read_aln.tags:
				if tag == 'X0': n += value
				elif tag == 'X1': n += value
				elif tag == 'XA': has_xa = True
			if n>1 and not has_xa:
				skipped_by_xa += 1
				# print("Skipping read %s/%s"%(read_aln.qname,'1' if read_aln.is_read1 else '2'),file=sys.stderr)
				continue
		if last_read_name == None:
			last_read_name = read_aln.qname
		if read_aln.qname != last_read_name:
			if options.verify_sortedness:
				if read_aln.qname in read_names:
					print("Error: Reads not grouped properly. Offending read:",read_aln.qname,file=sys.stderr)
					sys.exit(1)
				read_names.add(read_aln.qname)
			process_read_allpairs(last_read_name, inputfile.references, alignments1, alignments2, insert_length_distribution)
			alignments1, alignments2 = [], []
			alignments1_starts, alignments2_starts = set(), set()
			last_read_name = read_aln.qname
		def add_read(alignments, alignments_starts, read_aln):
			global skipped_duplicates
			if read_aln.pos in alignments_starts:
				skipped_duplicates += 1
				#print("Skipping",read_aln.qname, file=sys.stderr)
			else:
				alignments.append(read_aln)
				alignments_starts.add(read_aln.pos)
		if read_aln.is_read1:
			add_read(alignments1, alignments1_starts, read_aln)
		else:
			add_read(alignments2, alignments2_starts, read_aln)
	process_read_allpairs(last_read_name, inputfile.references, alignments1, alignments2, insert_length_distribution)
	if skipped_by_xa > 0:
		print("Skipped %d ambiguosly mapped reads for which no XA tag was present."%skipped_by_xa,file=sys.stderr)
	if skipped_duplicates > 0:
		print("Skipped %d duplicate alignments."%skipped_duplicates,file=sys.stderr)

if __name__ == '__main__':
	sys.exit(main())
