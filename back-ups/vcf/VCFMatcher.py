#!/usr/bin/env python
from __future__ import print_function, division
import vcf

from Variant.Indel import *

__author__ = "Louis Dijkstra"

"""
	VCFMatcher.py contains the class VCFMatcher that matches two VCF files together. 
"""


def returnAutosome (autosome):
	"""Returns an integer denoting the chromosome."""
	if len(autosome) > 3 and autosome[:3] == 'chr':
		autosome = autosome[3:]
	if autosome == 'X' or autosome == 'x':
		return 23
	if autosome == 'Y' or autosome == 'y':
		return 24	
	return int(autosome)

class VCFMatcher:
	
	def __init__(self, vcf_file, length_threshold = 100, centerpoint_dist_threshold = 100):
		self.vcf_reader 	= vcf.Reader(open(vcf_file))
		self.length_thres 	= length_threshold
		self.dist_thres 	= centerpoint_dist_threshold

	def retrieveNumberOfDeletions(self, list_of_autosomes = range(1,23)):
		"""Returns the number of deletions on the specified autosomes."""
		n_deletions = 0 		
		for vcf_record in self.vcf_reader:
			if not returnAutosome(vcf_record.CHROM) in list_of_autosomes:
				continue
			n_deletions += 1
		return n_deletions
			
	def retrieveMatchingDeletions(self, deletion, search_range = 10000):
		"""Returns a list of all matching VCF records in the VCF file."""

		# fetch all vcf_records from the same region
		autosome = 'chr' + str(returnAutosome(deletion.chromosome))
		start_search_range = max(0, deletion.start - search_range)
		end_search_range = deletion.end + search_range 
		candidate_vcf_records = self.vcf_reader.fetch(autosome, start = start_search_range, end = end_search_range)
		
		matching_deletions = [] 
		for candidate_vcf in candidate_vcf_records:
			candidate_deletion = Deletion(candidate_vcf)
			if deletion.similarTo(candidate_deletion, difference_length_threshold = self.length_thres, difference_centerpoint_threshold = self.dist_thres):
				matching_deletions.append(candidate_deletion)

		return matching_deletions 
	
		
