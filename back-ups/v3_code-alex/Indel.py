#!/usr/bin/env python
from __future__ import print_function, division
import vcf

__author__ = "Louis Dijkstra"

"""
	Indel.py contains the classes/functions for working with VCF files. 
	Used by the VAF Estimator and the Somatic Mutation Caller. 
"""

def isIndel (vcf_record):
	"""Determines whether a vcf record represents an indel or not. WARNING: only considers the first alternative (ALT[0]); others are neglected."""
	return (isDeletion(vcf_record) or isInsertion(vcf_record))

def isDeletion (vcf_record):
	"""Determines whether a vcf record represents a deletion or not. NOTE: only considers first alternative (ALT[0]); others are neglected"""
	if 'SVTYPE' in vcf_record.INFO: # checks whether the key SVTYPE is present in INFO
		if vcf_record.INFO['SVTYPE'][0] == 'DEL':	
			return True
	else: # makes statement on the basis of length of ALT and REF 
		if len(vcf_record.REF) > 1 and len(vcf_record.ALT[0]) == 1:
			return True
	return False

def isInsertion (vcf_record):
	"""Determines whether a vcf record represents an insertion or not. NOTE: only considers first alternative (ALT[0]); others are neglected"""
	if 'SVTYPE' in vcf_record.INFO: # checks whether the key SVTYPE is present in INFO
		if vcf_record.INFO['SVTYPE'][0] == 'INS':	
			return True
	else: # makes statement on the basis of length of ALT and REF 
		if len(vcf_record.REF) == 1 and len(vcf_record.ALT[0]) > 1:
			return True
	return False

def returnIndelLength (vcf_record):
	"""Returns variation length of an indel given its vcf record"""
	if 'SVLEN' in vcf_record.INFO:
		return abs(vcf_record.INFO['SVLEN'][0]) 
	else:
		return abs(len(vcf_record.REF) + 1 - len(vcf_record.ALT[0]))

def getDelta (vcf_record):
	"""Returns the delta (variant length) of an indel. Is positive in case of a deletion and negative for an insertion."""
	if isDeletion(vcf_record): 
		return returnIndelLength(vcf_record)
	return -1.0 * returnIndelLength(vcf_record)

class Indel: 
	"""Class for representing an indel (deletion/insertion)."""
	def __init__(self, vcf_record):
		self.vcf_record = vcf_record 
		self.chromosome = str(vcf_record.CHROM)
		self.length 	= returnIndelLength(vcf_record)

	def getTrueVAF (self):
		"""Returns true variant allele frequency as given in VCF file."""
		if self.vcf_record.heterozygosity == 0.5:
			return 0.5 
		else:
			return 1.0		

class Deletion(Indel):
	"""Class for representing a deletion. """
	def __init__(self, vcf_record):
		Indel.__init__(self, vcf_record)
		# Interval which is deleted [start, end]:
		self.start = self.vcf_record.POS + 1
		self.end   = self.start + self.length - 1
		# centerpoints are the points in the middle of the deletion; when length is even, there are two, otherwise one
		if self.length % 2 == 0: 
			self.centerpoints = [self.start + self.length / 2 - 1, self.start + self.length / 2]
		else: 
			self.centerpoints = [self.start + self.length / 2]

	def print(self):
		"""Prints most relevant data on the deletion to standard output"""
		print(self.chromosome, '\t', self.vcf_record.POS, '\t', self.length, '\t', end = "")

class Insertion(Indel):
	"""Class for represenitng an insertion."""
	def __init__(self, vcf_record):
		Indel.__init__(self, vcf_record)
		self.position = vcf_record.POS 

	def print(self):
		"""Prints most relevant data on the deletion to standard output"""
		print(self.chromosome, '\t', self.vcf_record.POS, '\t', self.length, '\t', end = "")
	
