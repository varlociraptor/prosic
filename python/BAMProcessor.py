#!/usr/bin/env python
from __future__ import print_function, division
import vcf
from collections import defaultdict
import pysam 
import numpy as np 

from Indel import *
from Alignments import *

__author__ = "Louis Dijkstra"


"""
	BAMProcessor.py contains the BAMProcessor superclass. It contains all the 
	functionality that all BAM Processors have in common, i.e., every other 
	BAMProcessor class is derived from this one. 
"""

class BAMProcessor: 
	"""Superclass for processing BAM files."""

	def __init__(self, bam_filename, search_range = 5000, primary_alignments_only = False):
		self.bam_reader 		= pysam.Samfile(bam_filename, "rb")
		self.search_range 		= search_range # range in which one searches for alignments
		self.primary_alignments_only 	= primary_alignments_only # when True, only primary alignments are considered
		
	def determineSupportDeletionSingleAlignment(self, deletion, alignment): 
		"Needs to be implemented in the subclass"
		pass 

	def determineSupportInsertionSingleAlignment(self, insertion, alignment): 
		"Needs to be implemented in the subclass"
		pass 


	def close(self): 
		"""Closes the BAM file."""
		self.bam_reader.close()

	def process(self, vcf_record):
		"""Collects the evidence (both overlapping and internal segment based) for a given indel (vcf record)."""
		if isDeletion(vcf_record):
			return processDeletion (self, Deletion(vcf_record))
		elif isInsertion(vcf_record):
			return processInsertion (self, Insertion(vcf_record))
		else:
			assert False 

	def processDeletion(self, deletion):
		"""Collects the evidence (both overlapping and internal segment based) for a given deletion."""
		isize, isize_prob, splits, splits_prob = [], [], [], []
		
		alignment_dict = defaultdict(list)

		# fetch the alignments in the vicinity of the deletion
		for alignment in self.bam_reader.fetch(deletion.chromosome, max(0, deletion.start - 1 - self.search_range), deletion.end + 1 + self.search_range):
			if alignment.isize == 0: # alignment is unmapped
				continue
			# determine whether this individual read is relevant 			
			if alignment.positions[0] < min(deletion.centerpoints) and alignment.positions[-1] > max(deletion.centerpoints): # read overlaps the centerpoints of the deletion 
				splits.append(self.determineSupportDeletionSingleAlignment(deletion, alignment))
				splits_prob.append(1.0 - convertPhredScore(alignment.mapq))
			else:
				alignment_dict[alignment.qname].append(alignment)

		# walk through the paired-end reads
		for qname, alignments in alignment_dict.iteritems():
			if len(alignments) == 2: # paired-end read 
				align_l, align_r = alignments[0], alignments[1]
				if align_r.positions[-1] < align_l.positions[0]:
					temp = align_r
					align_r = align_l
					align_l = temp
				interval_segm_start = align_l.positions[-1] + 1
				interval_segm_end  = align_r.positions[0] - 1
				if interval_segm_start <= min(deletion.centerpoints) and interval_segm_end >= max(deletion.centerpoints):
					isize.append(interval_segm_end - interval_segm_start + 1)
					isize_prob.append((1.0 - convertPhredScore(align_l.mapq)) * (1.0 - convertPhredScore(align_r.mapq)))

		return np.array(isize), np.array(isize_prob), np.array(splits), np.array(splits_prob)

	def processInsertion(self, insertion):
		"""Collects the evidence (both overlapping and internal segment based) for a given deletion."""
		isize, isize_prob, splits, splits_prob = [], [], [], []
		
		alignment_dict = defaultdict(list)

		# fetch the alignments in the vicinity of the deletion
		for alignment in self.bam_reader.fetch(insertion.chromosome, max(0, insertion.position - self.search_range), insertion.position + 1 + self.search_range):
			if alignment.isize == 0: # alignment is unmapped
				continue
			# determine whether this individual read is relevant 			
			if alignment.positions[0] < insertion.position and alignment.positions[-1] > insertion.position:
				splits.append(self.determineSupportInsertionSingleAlignment(insertion, alignment))
				splits_prob.append(1.0 - convertPhredScore(alignment.mapq))
			else:
				alignment_dict[alignment.qname].append(alignment)

		# walk through the paired-end reads
		for qname, alignments in alignment_dict.iteritems():
			if len(alignments) == 2: # paired-end read 
				align_l, align_r = alignments[0], alignments[1]
				if align_r.positions[-1] < align_l.positions[0]:
					temp = align_r
					align_r = align_l
					align_l = temp
				interval_segm_start = align_l.positions[-1] + 1
				interval_segm_end  = align_r.positions[0] - 1
				if interval_segm_start <= insertion.position and interval_segm_end >= insertion.position + 1:
					isize.append(interval_segm_end - interval_segm_start + 1)
					isize_prob.append((1.0 - convertPhredScore(align_l.mapq)) * (1.0 - convertPhredScore(align_r.mapq)))

		return np.array(isize), np.array(isize_prob), np.array(splits), np.array(splits_prob)

