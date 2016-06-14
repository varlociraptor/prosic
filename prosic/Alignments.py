#!/usr/bin/env python
from __future__ import print_function, division
import vcf
import pysam 

# import from Indel.py (the classes 'Indel', 'Deletion' and 'Insertion' are needed)
from Indel import *

__author__ = "Louis Dijkstra"

"""
	Alignments.py contains the classes for representing different types of alignments.  
	Used by the VAF Estimator and the Somatic Mutation Caller. 
"""

def convertPhredScore (phred_score):
	"""Converts phred score to probability"""
	return 10.0 ** (-float(phred_score) / float(10.0))

def returnMinimumDifference (list1, list2):
	"""Returns minimal difference between all possible pairs of elements contained in list1 and list2."""
	min_diff = float("inf")
	for e1 in list1:
		for e2 in list2:
			if min_diff > abs(e1 - e2):
				min_diff = abs(e1 - e2)
	return min_diff

class DummyAlignment: 
	"""Class that only contains the observation (value) and alignment probability."""
	def __init__(self, value, probability):
		self.value = value
		self.probability = probability
	
class PairedEndAlignment:
	"""Class for alignments of paired-end reads."""
	def __init__(self, alignment1, alignment2):
		# determine which alignments lies left of the other
		if alignment1.positions[-1] < alignment2.positions[0]:
			self.alignment_left = alignment1 
			self.alignment_right = alignment2 
		else:
			self.alignment_left = alignment2 
			self.alignment_right = alignment1
		
		self.internal_segment_start = self.alignment_left.positions[-1] + 1
		self.internal_segment_end = self.alignment_right.positions[0] - 1 
		
		self.value = self.alignment_right.positions[0] - self.alignment_left.positions[-1] + 1 # internal segment length
		self.probability = self.returnAlignmentProbability()

	def returnAlignmentProbability(self): 
		"""Returns alignment probability of paired-end read. Assumes independence"""
		return (1.0 - convertPhredScore(self.alignment_left.mapq)) * (1.0 - convertPhredScore(self.alignment_right.mapq))

	def relevantForDeletion (self, deletion):
		"""A paired end alignment is deemed relevant for a deletion when the left alignment lies left of 
		   the centerpoint of the deletion, while the right alignment lies right."""
		if self.internal_segment_start <= min(deletion.centerpoints) and self.internal_segment_end >= max(deletion.centerpoints):
			return True
		return False 
		
	def relevantForInsertion (self, insertion): 
		"""A paired end alignment is deemed relevant for an insertion when the left alignment lies left of 
		   the breakpoint, while the right alignment lies right."""
		if self.internal_segment_start <= insertion.position and self.internal_segment_end >= insertion.position + 1: 
			return True
		return False
	
	def print(self):
		"""Prints the paired-end alignment to std out."""
		print('Paired-end alignment |\tLeft: [', self.alignment_left.positions[0], ',', self.alignment_left.positions[-1], ']\tRight: [', self.alignment_right.positions[0], ',', self.alignment_right.positions[-1], ']\t| Internal segm. length: ', self.value, '\t| Prob.: ', self.probability) 
		
		
class OverlappingAlignment: 
	"""Class for an alignment of one read (no-paired end)"""
	def __init__(self, alignment): 
		self.alignment 		= alignment 
		self.start 		= alignment.positions[0]
		self.end 		= alignment.positions[-1]
		self.probability 	= (1.0 - convertPhredScore(self.alignment.mapq)) 
		self.value 		= None 	# 1 when supports presence indel, 0 otherwise
		self.splits 		= self.returnSplits()

	def returnSplits(self):
		"""Returns a list of the splits present in the alignment."""
		splits = [] 

		i = self.start + 1
		while i <= self.end:
			if i in self.alignment.positions: # position is present in alignment -> no split
				i = i + 1
			else:
				# start of split
				start_split = i 
				i = i + 1
				while i not in self.alignment.positions:
					i = i + 1
				end_split = i - 1
				i = i + 1
				splits.append([start_split, end_split])
		return splits 

	def relevantForDeletion(self, deletion):
		"""Alignment is considered relevant for a given deletion when the alignment overlaps the centerpoint of the deletion.""" 
		if self.start > max(deletion.centerpoints) or self.end < min(deletion.centerpoints):
			return False 
		self.determineSupportDeletion(deletion)
		return True 
		
	def determineSupportDeletion(self, deletion):
		"""Determines whether the alignment supports the presence/absence of a given deletion."""
		self.value = 0 # initially 0. Supports deletion when proof is there 

		# walk through all splits and see whether there is one that supports the presence of the deletion
		for split in self.splits: 
			length_split = split[1] - split[0] + 1
			if abs(length_split - deletion.length) > 20:
				continue
			centerpoints_split = None 
			if length_split % 2 == 0: # when length is even, there are two centerpoints
				centerpoints_split = [split[0] + length_split / 2 - 1, split[0] + length_split / 2]
			else:
				centerpoints_split = [split[0] + length_split / 2]
			if returnMinimumDifference(centerpoints_split, deletion.centerpoints) <= 50:
				self.value = 1 
		
	def relevantForInsertion (self, insertion):
		"""Alignment is considered relevant for a given insertion when the alignment overlaps the position where the insertion occurred."""
		if self.start > insertion.position or self.end < insertion.position + 1:
			return False
		self.determineSupportInsertion(insertion)
		return True

	def determineSupportInsertion(self, insertion):
		"""Determines whether the alignment supports the presence/absence of a given insertion."""
		self.value = 0 # initially 0. Support insertion only when proof is there

		cigar_line = self.alignment.cigar 
		
		# walk through the read and find a insertion splits
		i = self.start
		for (cigar_type, cigar_length) in cigar_line:
			if cigar_type == 1: # insertion! 
				if abs(cigar_length - insertion.length) > 20:
					continue
				centerpoints_split = None
				if cigar_length % 2 == 0: 
					centerpoints_split = [i + cigar_length / 2 - 1, i + cigar_length / 2]
				else:
					centerpoints_split = [i + cigar_length / 2]
				if returnMinimumDifference(centerpoints_split, [insertion.pos, insertion.pos + 1]) <= 50:
					self.value = 1
	
