#!/usr/bin/env python
from __future__ import print_function, division
import vcf
from collections import defaultdict
import pysam 

from Indel import *
from Alignments import *

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
		self.paired_end_alignments 	= [] # paired-end alignments retrieved
		self.overlapping_alignments 	= [] # overlapping alignments retrieved
		
	def reset(self): 
		"""Resets the results found for the previous search."""
		self.paired_end_alignments = [] 
		self.overlapping_alignments = []	
		
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

	def fetchAlignments (self, chromosome, start, end): 
		"""Returns all alignments from a the range [start - search_range, end + search_range] in a default dictionary, where keys are the qnames of the alignments"""
		alignment_dict = defaultdict(list)
		# check whether chromosome has the right structure
		chromosome = self.getChromosome(chromosome)
		if self.primary_alignments_only:
			for alignment in self.bam_reader.fetch(self.getChromosome(chromosome), max(0, start - self.search_range), end + self.search_range):
				if alignment.is_unmapped or alignment.is_secondary: # aligned is not mapped or not primary
					continue
				alignment_dict[alignment.qname].append(alignment)
		else:
			for alignment in self.bam_reader.fetch(self.getChromosome(chromosome), max(0, start - self.search_range), end + self.search_range):
				if alignment.is_unmapped: # aligned is not mapped or not primary
					continue
				alignment_dict[alignment.qname].append(alignment)
		return alignment_dict

	def process(self, vcf_record):
		"""Collects the evidence (both overlapping and internal segment based) for a given indel (vcf record)."""
		if isDeletion(vcf_record):
			return processDeletion (self, Deletion(vcf_record))
		elif isInsertion(vcf_record):
			return processInsertion (self, Insertion(vcf_record))
			
	def processDeletion (self, deletion): 
		"""Collects the evidence (both overlapping and internal segment based) for a given deletion."""
		self.reset()
		
		# obtain all the potentially relevant alignments
		alignment_dict = self.fetchAlignments(deletion.chromosome, deletion.start - 1, deletion.end + 1)
		
		for qname, alignments in alignment_dict.iteritems(): # walk through all grouped alignments
			# test whether they are applicable for the given deletion
			if len(alignments) == 2: # paired-end read 
				paired_end_alignment = PairedEndAlignment(alignments[0], alignments[1])
				if paired_end_alignment.relevantForDeletion(deletion): 
					self.paired_end_alignments.append(paired_end_alignment)
				else: # potential overlapping alignments
					overlapping_alignment1 = OverlappingAlignment(alignments[0])
					overlapping_alignment2 = OverlappingAlignment(alignments[1])
					if overlapping_alignment1.relevantForDeletion(deletion):
						self.overlapping_alignments.append(overlapping_alignment1)
					if overlapping_alignment2.relevantForDeletion(deletion):
						self.overlapping_alignments.append(overlapping_alignment2)
			elif len(alignments) == 1: # potential overlapping alignment
				overlapping_alignment = OverlappingAlignment(alignments[0])
				if overlapping_alignment.relevantForDeletion(deletion):
					self.overlapping_alignments.append(overlapping_alignment)
			
		return self.paired_end_alignments, self.overlapping_alignments	
			
	def processInsertion (self, insertion): 
		"""Collects the evidence (both overlapping and internal segment based) for a given insertion."""
		self.reset()
		
		# obtain all the potentially relevant alignments		
		alignment_dict = self.fetchAlignments(insertion.chromosome, insertion.position, insertion.position + 1)
		
		for qname, alignments in alignment_dict.iteritems(): # walk through all grouped alignments
			# test whether they are applicable for the given insertion
			
			if len(alignments) == 2: # paired-end read 
				paired_end_alignment = PairedEndAlignment(alignments[0], alignments[1])
				if paired_end_alignment.relevantForInsertion(insertion): 
					self.paired_end_alignments.append(paired_end_alignment)
				else: 
					overlapping_alignment1 = OverlappingAlignment(alignments[0])
					overlapping_alignment2 = OverlappingAlignment(alignments[1])
					if overlapping_alignment1.relevantForInsertion(insertion):
						self.overlapping_alignments.append(overlapping_alignment1)
					if overlapping_alignment2.relevantForInsertion(insertion):
						self.overlapping_alignments.append(overlapping_alignment2)
			elif len(alignments) == 1: # potential overlapping alignment
				overlapping_alignment = OverlappingAlignment(alignments[0])
				if overlapping_alignment.relevantForInsertion(insertion):
					self.overlapping_alignments.append(overlapping_alignment)
			
		return self.paired_end_alignments, self.overlapping_alignments
		
	

