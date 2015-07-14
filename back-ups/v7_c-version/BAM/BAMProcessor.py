#!/usr/bin/env python
from __future__ import print_function, division
import vcf
from collections import defaultdict
import pysam 
import numpy as np 

from Variant.Indel import *
from BAM.Alignments import *

__author__ = "Louis Dijkstra"

"""
	BAMProcessor.py contains the functionality to process a given BAM file. 
"""

class BAMProcessor: 
	"""Class for processing BAM files. Contains all functionality needed for retrieving alignments."""
	def __init__(self, bam_filename, search_range = 5000, primary_alignments_only = False):
		self.bam_reader 		= pysam.Samfile(bam_filename, "rb")
		self.search_range 		= search_range # range in which one searches for alignments
		self.primary_alignments_only 	= primary_alignments_only # when True, only primary alignments are considered
			
	def close(self): 
		"""Closes the BAM file."""
		self.bam_reader.close()

	def getChromosome(self, chromosome):
		"""Returns the chromosome description needed by pysam"""
		if len(chromosome) < 3: 
			return 'chr' + chromosome
		if chromosome[:3] == 'chr':
			return chromosome
		return 'chr' + chromosome

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
		for alignment in self.bam_reader.fetch(self.getChromosome(deletion.chromosome), max(0, deletion.start - 1 - self.search_range), deletion.end + 1 + self.search_range):
			if alignment.isize == 0: # alignment is unmapped
				continue
			# determine whether this individual read is relevant 			
			if alignment.positions[0] < min(deletion.centerpoints) and alignment.positions[-1] > max(deletion.centerpoints):
				# add alignment observation 
				observation = 0 # initially 0. Supports deletion when proof is there
				i = alignment.positions[0] + 1
				while i <= alignment.positions[-1]:
					if i in alignment.positions: # position is present in alignment -> no split
						i += 1
					else:
						start_split = i
						i += 1
						while i not in alignment.positions:
							i += 1
						end_split = i-1
						i += 1
						length_split = end_split - start_split + 1
						if abs(length_split + 1 - deletion.length) > 20:
							continue
						centerpoints_split = None 
						if length_split % 2 == 0: # when length is even, there are two centerpoints
							centerpoints_split = [start_split + length_split / 2 - 1, start_split + length_split / 2]
						else:
							centerpoints_split = [start_split + length_split / 2]
						if returnMinimumDifference(centerpoints_split, deletion.centerpoints) <= 50:
							observation = 1
				splits.append(observation)
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
		for alignment in self.bam_reader.fetch(self.getChromosome(insertion.chromosome), max(0, insertion.position - self.search_range), insertion.position + 1 + self.search_range):
			if alignment.isize == 0: # alignment is unmapped
				continue
			# determine whether this individual read is relevant 			
			if alignment.positions[0] < insertion.position and alignment.positions[-1] > insertion.position:
				# add alignment observation 
				observation = 0 # initially 0. Supports insertion when proof is there
				cigar_line = alignment.cigar 
				# walk through the read and find a insertion splits
				i = alignment.positions[0]
				for (cigar_type, cigar_length) in cigar_line:
					if cigar_type == 1: # insertion! 
						if abs(cigar_length - insertion.length) > 20:
							continue
						centerpoints_split = None
						if cigar_length % 2 == 0: 
							centerpoints_split = [i + cigar_length / 2 - 1, i + cigar_length / 2]
						else:
							centerpoints_split = [i + cigar_length / 2]
						if returnMinimumDifference(centerpoints_split, [insertion.position, insertion.position + 1]) <= 50:
							observation = 1
				splits.append(observation)
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

