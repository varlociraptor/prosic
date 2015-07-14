#!/usr/bin/env python
from __future__ import print_function, division
from BAMProcessor import * # import the super class 

__author__ = "Louis Dijkstra"

"""
	BWABAMProcessor.py contains the functionality to process a given BWA(-memm) BAM file. 
"""

class BWABAMProcessor(BAMProcessor): 
	"""Class for processing BWA BAM files. Contains all functionality needed for retrieving alignments."""
	def __init__(self, bam_filename, search_range = 5000, primary_alignments_only = False):
		BAMProcessor.__init__(self, bam_filename, search_range = search_range, primary_alignments_only = primary_alignments_only)	

	def determineSupportDeletionSingleAlignment(self, deletion, alignment): 
		for (cigar_type, cigar_length) in alignment.cigar: # walk through the cigar string
			if cigar_type == 2: # deletion
				if abs(cigar_length - deletion.length) <= 1: 
					return 1 
		return 0 

	def determineSupportInsertionSingleAlignment(self, insertion, alignment):
		for (cigar_type, cigar_length) in alignment.cigar: # walk through the cigar string
			if cigar_type == 1: # insertion
				if abs(cigar_length - insertion.length) <= 1: 
					return 1 
		return 0 
	

