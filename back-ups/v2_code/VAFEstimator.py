#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
from collections import defaultdict
from scipy.stats import norm 
from scipy.optimize import fminbound
from scipy.optimize import bisect
import os
import sys
import vcf
import pysam
import math

__author__ = "Louis Dijkstra"

usage = """%prog [options] <mu> <sigma> <vcf-filename> <bam-filename> 

	<mu>		the mean of the null normal 
				distribution for the internal 
				segment length when not 
				affected by an indel. 
	<sigma>		the standard deviation of the 
				null normal distribution for 
				the internal segment length 
				when not affected by an indel
	<vcf-file> 	tabix-indexed VCF file containing the 
				indels which allele frequency 
				needs to be estimated. 
				NOTE: File needs to be sorted 
				and indexed. 
	<bam-file>	BAM file containing the 
				alignments. 
				NOTE:File needs to be sorted 
				and indexed. 

NOTE: the current implementation
		1) 	ignores indels in the vcf file 
			for which more than one 
			alternative (ALT) are provided. 	
		2)	only considers primary 
			alignments; secondary/ternary etc. 
			alignments are neglected. 
"""


def convertPhredScore (phred_score):
	"""Converts phred score to probability"""
	return 10.0 ** (-phred_score / 10.0)
	
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

def isIndel (vcf_record):
	"""Determines whether a vcf record represents an indel or not. WARNING: only considers the first alternative (ALT[0]); others are neglected."""
	return (isDeletion(vcf_record) or isInsertion(vcf_record))

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
	
def returnEpsilon0(delta):
	return 0.0

def returnEpsilon1(delta):
	return 0.0


def returnChromosome(chromosome):
	chromosome = str(chromosome)
	if chromosome == 'X' or chromosome == 'Y':
		return chromosome
	#if chromosome[:3] == 'chr':
	#	return chromosome
	#return "chr" + chromosome

def returnDistance (list1, list2): 
	"""
		Returns minimal difference between all possible pairs of elements contained
		in two different lists 
	"""
	min_dist = float("inf")
	for i in range(len(list1)):
		for j in range(len(list2)):
			if abs(list1[i] - list2[j]) < min_dist:
				min_dist = abs(list1[i] - list2[j])
	
	return min_dist

class Indel: 
	"""
		Class for representing an indel (deletion/insertion). 
	"""
	def __init__(self, vcf_record):
		self.vcf_record = vcf_record 
		self.chromosome = returnChromosome(vcf_record.CHROM)
		self.length 	= returnIndelLength(vcf_record)

	def getTrueVAF (self):
		"""Returns true variant allele frequency as given in VCF file"""
		if self.vcf_record.heterozygosity == 0.5:
			return 0.5 
		else:
			return 1.0		

class Deletion(Indel):
	"""
		Class for representing a deletion. 
	"""
	def __init__(self, vcf_record):
		Indel.__init__(self, vcf_record)
		self.start = self.vcf_record.POS + 1
		self.end = self.start + self.length
		if self.length % 2 == 0: # when length is even, there are two 'centerpoints'
			self.centerpoints = [self.start + self.length / 2 - 1, self.start + self.length / 2]
		else: 
			self.centerpoints = [self.start + self.length / 2]

	def print(self):
		print('Deletion  | chr: ', self.chromosome, ' |\tfrom: ', self.start, '\tto: ', self.end, '\t | length: ', self.length) 

class Insertion(Indel):
	"""
		Class for represenitng an insertion.
	"""
	def __init__(self, vcf_record):
		Indel.__init__(self, vcf_record)
		self.position = vcf_record.POS 

	def print(self):
		print('Insertion | chr: ', self.chromosome, ' |\tposition: ', self.position, '\t\t\t | length: ', self.length) 
	

class PairedEndAlignment:
	"""
		Class for paired end alignments. Contains two alignments of a paired-end read. 
	"""
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
		self.determineAlignmentProbability()


	def determineAlignmentProbability(self): 
		"""
			Returns alignment probability of paired-end read.
		"""
		self.probability =  (1.0 - convertPhredScore(self.alignment_left.mapq)) * (1.0 - convertPhredScore(self.alignment_right.mapq))


	def relevantForDeletion (self, deletion):
		"""
			A paired end alignment is relevant for a deletion when the left alignment lies left of 
			the centerpoint of the deletion, while the right alignment lies right. 
		"""
		if len(deletion.centerpoints) == 2: 
			if self.internal_segment_start <= deletion.centerpoints[0] and self.internal_segment_end >= deletion.centerpoints[1]:
				return True
		else: 
			if self.internal_segment_start <= deletion.centerpoints[0] and self.internal_segment_end >= deletion.centerpoints[0]:
				return True
		return False	
		
	def relevantForInsertion (self, insertion):
		if self.internal_segment_start <= insertion.position and self.internal_segment_end >= insertion.position + 1: 
			return True
		return False
	
	def print(self):
		print('Paired-end alignment |\tLeft: [', self.alignment_left.positions[0], ',', self.alignment_left.positions[-1], ']\tRight: [', self.alignment_right.positions[0], ',', self.alignment_right.positions[-1], ']\t| Internal segm. length: ', self.value, '\t| Prob.: ', self.probability) 
		
		
class OverlappingAlignment: 
	"""
		Class for one alignment.
	"""
	def __init__(self, alignment): 
		self.alignment 	= alignment 
		self.start 	= alignment.positions[0]
		self.end 	= alignment.positions[-1]
		self.determineAlignmentProbability()
		self.value 	= None 	# 1 when supports presence indel, 0 otherwise

	def determineAlignmentProbability(self):
		"""
			Determines alignment probability of this overlapping read.
		"""
		self.probability = (1.0 - convertPhredScore(self.alignment.mapq)) 



	def determineSplits (self):
		"""
			Determines the splits in the alignments.
		"""
		self.splits = [] 

		i = self.start + 1

		while i <= self.end:
			if i in self.alignment.positions: # position is present in alignment -> no split
				i = i + 1
			else:
				# start of split!
				start_split = i 
				i = i + 1
				while i not in self.alignment.positions:
					i = i + 1
				end_split = i - 1
				i = i + 1
				self.splits.append([start_split, end_split])
	
			

	def relevantForDeletion (self, deletion):
		"""Determines whether the alignment provides evidence for or against a given deletion"""
		if len(deletion.centerpoints) == 2: 
			if self.start > deletion.centerpoints[1] or self.end < deletion.centerpoints[0]:
				return False
		else: 
			if self.start > deletion.centerpoints[0] or self.end < deletion.centerpoints[0]:
				return False
		
		self.determineSupportDeletion(deletion)
		return True 
		
	def determineSupportDeletion(self,deletion):
		"""Determines the evidence for or against a given deletion"""

		self.value = 0 # initially 0. Supports deletion when proof is there 

		# determine all the splits 
		self.determineSplits () 

		# walk through all splits and see whether there is one that supports the presence of the deletion
		for split in self.splits: 
			start_split = split[0]
			end_split = split[1]
			length_split = end_split - start_split + 1
			
			if abs(length_split - deletion.length) <= 20:
				 # determine center points split 
				centerpoint_split = None
				if length_split % 2 == 0: # when lenght is even, there are two 'centerpoints'
					centerpoint_split = [start_split + length_split / 2 - 1, start_split + length_split / 2]
				else: 
					centerpoint_split = [start_split + length_split / 2]

				if returnDistance(centerpoint_split, deletion.centerpoints) <= 50: 
					self.value = 1 
				else:
					self.probability = 0.0 
			else:
				self.probability = 0.0 

	
		
	def relevantForInsertion (self, insertion):
		if self.start > insertion.position or self.end <= insertion.position: # alignment does not overlap
			return False 
		if isSplitRead(self.alignment):
			self.value = 1 
		else:
			self.value = 0 
		return True



class BAMProcessor: 
	"""
		Class for processing BAM files. Contains all functionality
		needed for retrieving alignments 
	"""
	
	def __init__(self, bam_filename, search_range=5000):
		self.bam_reader = pysam.Samfile(bam_filename, "rb")
		self.search_range = search_range
		self.paired_end_alignments = [] 
		self.overlapping_alignments = [] 
		
	def reset(self): 
		self.paired_end_alignments = [] 
		self.overlapping_alignments = []	
		
	def close(self): 
		self.bam_reader.close()
	
	def fetchAlignments (self, chromosome, start, end): 
		"""Returns all alignments from a the range [start - search_range, end + search_range] in a default dictionary, where keys are the qnames of the alignments"""
		alignment_dict = defaultdict(list)
		for alignment in self.bam_reader.fetch(chromosome, max(0, start - self.search_range), end + self.search_range):
			if not alignment.is_unmapped: # alignment is mapped
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
				else: 
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
		
	



class TestCase: 
	"""
		Class to store all results for one test case (deletion/insertion)
	"""
	def __init__(self, is_deletion, mle_disc, unique_disc, mle_cont, unique_cont, ci, indel, paired_end_alignments, overlapping_alignments, posterior_distr):
		self.is_deletion 		= is_deletion		
		self.mle_disc 			= mle_disc
		self.unique_disc 		= unique_disc
		self.mle_cont 			= mle_cont
		self.unique_cont 		= unique_cont
		self.ci 			= ci 
		self.true_vaf 			= indel.getTrueVAF() 		
		self.position 			= indel.vcf_record.POS 
		self.length 			= indel.length
		self.n_paired_end_alignments 	= len(paired_end_alignments)
		self.n_overlapping_alignments 	= len(overlapping_alignments)
		self.posterior_distr 		= posterior_distr 

		self.X = []
		self.Y = []

		for a in paired_end_alignments: 
			self.X.append(a.value)

		for a in overlapping_alignments:
			self.Y.append(a.value)

		if unique_disc:
			if mle_disc == self.true_vaf:
				self.correct_disc = True
			else:
				self.correct_disc = False

		# continuous estimate is considered correct when the 95% CI contains the true VAF 
		if ci[0] <= self.true_vaf and ci[1] >= self.true_vaf:
			self.correct_cont = True
		else:
			self.correct_cont = False

	def print(self, min_length = 0):
		"""
			Prints the relevant data of this test case to the command line
		"""
		if self.is_deletion and self.length >= min_length: # In case of a deletion
			print('--------------------------------------------------------------------------------------------------------')
			print('DEL\t|\t', self.true_vaf, '\t|\t', self.position, '\t', self.length) # print info on deletion 
			print('\nAvailable data: ')
			print('\t\t Internal segment lengths :', self.n_paired_end_alignments, '\t', self.X) 
			print('\t\t Overlapping alignments   :', self.n_overlapping_alignments, '\t', self.Y)
			print('\nDiscrete estimate: ')
			if self.unique_disc:
				if self.correct_disc:
					print('\t\tunique\tcorrect\t\t|\t', self.mle_disc, '\t|\t', max(self.posterior_distr), '\t', self.posterior_distr)
				else:
					print('\t\tunique\tincorrect\t|\t', self.mle_disc, '\t|\t', max(self.posterior_distr), '\t', self.posterior_distr)
			else:
				print('\t\tnot unique') 
			print('\nContinuous estimate: ')
			if self.unique_cont:
				if self.correct_cont:
					print('\t\tunique\tcorrect\t\t|\t', self.mle_cont, '\t|\t', self.ci)
				else:
					print('\t\tunique\tincorrect\t|\t', self.mle_cont, '\t|\t', self.ci)
			else:
				print('\t\tnot unique') 
			print('--------------------------------------------------------------------------------------------------------')
			print('\n\n')
					
		elif self.length >= min_length:
			if self.unique:
				if self.correct:
					print('INS\tunique\tV\t', self.position, '\t', self.mle, '\t', self.true_vaf,'\t', self.length, '\t', max(self.posterior_distr), '\t',self.posterior_distr, '\t', self.n_paired_end_alignments, '\t', self.n_overlapping_alignments, '\t', self.X, '\t', self.Y) 
				else:
					print('INS\tunique\tX\t', self.position, '\t', self.mle, '\t', self.true_vaf,'\t', self.length, '\t', max(self.posterior_distr), '\t', self.posterior_distr, '\t', self.n_paired_end_alignments, '\t', self.n_overlapping_alignments, '\t', self.X, '\t', self.Y) 		

			else:
				print('INS\tnot unique\t\t', self.position, '\t', self.mle, '\t', self.true_vaf,'\t', self.length, '\t',max(self.posterior_distr), '\t', self.posterior_distr, '\t', self.n_paired_end_alignments, '\t', self.n_overlapping_alignments) 
					


	def confidentEnough (self, threshold):
		"""
			Returns True when there is a posterior probability larger or equal to the given threshold; otherwise False. 
		"""
		for p in self.posterior_distr:
			if p >= threshold:
				return True
		return False

	def enoughEvidence (self, n):
		"""
			Returns true when the number of relevant reads is equal or exceeds a given n
		"""
		if self.n_paired_end_alignments + self.n_overlapping_alignments >= n:
			return True
		return False

	def isOfLength(self, a, b): 
		"""Determines whether this test case falls within a certain length class"""
		if a <= self.length <= b:
			return True
		return False


		

class TestDataCollector: 
	"""
		Only used for testing. TODO: Remove later
	"""
	def __init__(self):
		self.testcases = [] 
		self.length_range = [[10, 29],[30,49],[50,99],[100,249]]
		self.min_length = 10 # TODO change?
		
	def updateDeletion(self, MLE_disc, unique_disc, MLE_cont, unique_cont, CI, deletion, paired_end_alignments, overlapping_alignments, posterior_distr):	
		self.testcases.append(TestCase(True, MLE_disc, unique_disc, MLE_cont, unique_cont, CI, deletion, paired_end_alignments, overlapping_alignments, posterior_distr))
		self.testcases[-1].print(min_length = self.min_length)

	def updateInsertion(self, MLE_disc, unique_disc, MLE_cont, unique_cont, CI, insertion, paired_end_alignments, overlapping_alignments, posterior_distr):	
		self.testcases.append(TestCase(False, MLE_disc, unique_disc, MLE_cont, unique_cont, CI, insertion, paired_end_alignments, overlapping_alignments, posterior_distr))
		#self.testcases[-1].print(min_length = self.min_length)


	def printTable3x3 (self, table, title):
		row_total = [0,0,0]
		column_total = [0,0,0]
		n = 0 
		# compute marginals
		for r in [0,1,2]:
			for c in [0,1,2]:
				row_total[r] = row_total[r] + table[r][c]
				column_total[c] = column_total[c] + table[r][c]
				n = n + table[r][c]

		correct = 0 
		for i in [0,1,2]:
			correct = correct + table[i][i]

		incorrect = n - correct

		print('------------------------TABLE', title, '------------------------')
		print('\t\t\t\t\t\ttruth')
		print('\t\t\t|\t0.0\t|\t0.5\t|\t1.0\t|\ttotal')
		print('\t--------------------------------------------------------------------------------------------')
		print('\t\t0.0\t|\t',table[0][0], '\t|\t', table[0][1],'\t|\t', table[0][2],'\t|\t', row_total[0])
		print('assign.\t\t0.5\t|\t',table[1][0], '\t|\t', table[1][1],'\t|\t', table[1][2],'\t|\t', row_total[1])
		print('\t\t1.0\t|\t',table[2][0], '\t|\t', table[2][1],'\t|\t', table[2][2],'\t|\t', row_total[2])
		print('\t--------------------------------------------------------------------------------------------')
		print('\t\ttotal\t|\t', column_total[0], '\t|\t', column_total[1], '\t|\t', column_total[2], '\t|\t',  n)
		print('\nCORRECT:\t', correct, ' /', correct / float(n) * 100.0, '%')
		print('INCORRECT:\t', incorrect, ' /', incorrect / float(n) * 100.0, '%\n')

	def summarize(self, confidence_level = 0.0):
		"""Prints a summary of the test results to the command line"""
		for range in self.length_range:
			cont_correct_case = 0.0 
			total = 0.0 
			
			table = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
			for test_case in self.testcases:
				if test_case.isOfLength(range[0], range[1]):
					total = total + 1
					if test_case.unique_disc and test_case.isOfLength(range[0], range[1]) and test_case.confidentEnough(confidence_level):
						i = int(test_case.mle_disc * 2.0)
						j = int(test_case.true_vaf * 2.0)
						table[i][j] = table[i][j] + 1
					if test_case.correct_cont:
						cont_correct_case = cont_correct_case + 1
			title = 'results for indels of lengths between ' + str(range[0]) + ' and ' + str(range[1])
			self.printTable3x3(table, title)
			print('\nContinuous case: in', int(cont_correct_case), 'of', int(total), 'cases, the 95% confidence interval contained the true VAF, i.e.,', cont_correct_case / total * 100, '%\n') 


				
class VAFLikelihoodMaximizer: 
	
	def __init__(self, mu, sigma, vaf_range = [0.0, 0.5, 1.0]):
		self.mu 	= float(mu) 
		self.sigma 	= float(sigma)
		self.vaf_range = vaf_range
		self.paired_end_alignments = []
		self.overlapping_alignments = []
		self.delta = 0.0
		self.epsilon0 = 0.0 
		self.epsilon1 = 0.0
	
	def setVAFRange (self, vaf_range):
		"""
			Sets the VAF range. If vaf_range = [] or 'None', the vaf_range is considered lie within 
			the unit interval. The MLE can, thus, range continuosly over that region.
		"""
		self.vaf_range = vaf_range

	def print(self): 
		"""Outputs info on the likelihood maximizer to the command line"""

		print('-------------------VAF Likelihood Maximizer-------------------')
		print('Mean (mu)                 :', self.mu)	
		print('Standard deviation (sigma):', self.sigma)
		if self.vaf_range == [] or self.vaf_range == None: 
			print('VAF support               : continuous (unit interval)') 
		else:	
			print('VAF support               :', self.vaf_range)
		print('--------------------------------------------------------------')
		
		
	def f(self, x, mean, std):
		"""
			Returns probability of observing x given a discretized normal distribution with mean mu and std of
			sigma, truncated below 0
		"""
		return (norm.cdf((x + 1.0 - mean) / std) - norm.cdf((x - mean) / std)) / (1.0 - norm.cdf(-1 * mean / std))
	
	def loglikelihood (self, vaf):
		logl = 0.0 
		for paired_end_alignment in self.paired_end_alignments: # internal segment based evidence
			x = paired_end_alignment.value
			p = paired_end_alignment.probability
			logl = logl + math.log(p * (vaf * self.f(x, self.mu + self.delta, self.sigma) + (1 - vaf) * self.f(x, self.mu, self.sigma)) + (1.0 - p))
		for overlapping_alignment in self.overlapping_alignments: # overlapping read evidence
			y = overlapping_alignment.value
			p = overlapping_alignment.probability 
			if y == 1: 
				logl = logl + math.log(p * (vaf * (1.0 - self.epsilon1) + (1.0 - vaf) * self.epsilon0) + (1.0 - p))
			else: 
				logl = logl + math.log(p * (vaf * self.epsilon1 + (1.0 - vaf) * (1.0 - self.epsilon0)) + (1.0 - p))
		return logl 

	def set(self, paired_end_alignments, overlapping_alignments, delta):
		self.paired_end_alignments = paired_end_alignments
		self.overlapping_alignments = overlapping_alignments
		self.delta = float(delta)
		self.epsilon0 = returnEpsilon0(delta)
		self.epsilon1 = returnEpsilon1(delta)

	def maximize (self, vcf_record, paired_end_alignments, overlapping_alignments):
		delta = getDelta (vcf_record)
		
		if self.vaf_range == None or self.vaf_range == []: # VAF lies within [0,1]
			return self.computeMLE_continuous (paired_end_alignments, overlapping_alignments, delta)  
		else:
			return self.computeMLE_discrete (paired_end_alignments, overlapping_alignments, delta)
			
	def computeMLE_discrete(self, paired_end_alignments, overlapping_alignments, delta): 
		"""
			Computes the MLE when the support of the VAF is discrete. 
			Returns both the MLE and a boolean, which is true when the 
			MLE is unique and false otherwise.
		"""
		self.set(paired_end_alignments, overlapping_alignments, delta)
		
		MLE = [self.vaf_range[0]]
		max_loglikelihood = self.loglikelihood (MLE[0])
		
		for i in range(1, len(self.vaf_range)):
			vaf = self.vaf_range[i]
			logl = self.loglikelihood (vaf) 
			if logl > max_loglikelihood: 
				max_loglikelihood = logl 
				MLE = [vaf]
			elif logl == max_loglikelihood: 
				MLE.append(vaf)
		
		if len(MLE) == 1: 
			return MLE[0], True
		else: 
			return MLE, False
	
	def globalMaximumExists (self):
		"""Determines whether a global maximum exists (for the continuous case)"""
		for a in self.paired_end_alignments: 
			if a.probability > 0.0 and (self.f(a.value, self.mu + self.delta, self.sigma) is not self.f(a.value, self.mu, self.sigma)): 
				return True
		for a in self.overlapping_alignments: 
			if a.probability > 0.0: 
				return True
		return False
	
	def help_loglikelihood(self, vaf):
		return -1 * self.loglikelihood(vaf)
	
	def computeMLE_continuous (self, paired_end_alignments, overlapping_alignments, delta):
		"""
			Computes the MLE when the support of the VAF is the unit interval. 
			Returns both the MLE and a boolean, which is true when the 
			MLE is unique and false otherwise.
		"""
		self.set(paired_end_alignments, overlapping_alignments, delta)
		if self.globalMaximumExists() : 
			return fminbound(self.help_loglikelihood, 0.0, 1.0, [], xtol=0.0000001), True # TODO check! 
		else: 
			return None, False
		

	def help_ci_loglikelihood(self, vaf):
		return self.loglikelihood(vaf) - self.max_logl + 1.92


	def computeMLE_95CI (self, vcf_record, paired_end_alignments, overlapping_alignments):
		"""
			Computes the MLE when the support of the VAF is the unit interval. 
			In addition, the approx. 95% confidence interval is given. 
			Output:
				- Maximum likelihood estimate
				- Boolean; true, when MLE is unique and false otherwise
				- 95% confidence interval
		"""
		delta = getDelta (vcf_record)
		self.set(paired_end_alignments, overlapping_alignments, delta)

		CI = [0.0, 1.0] # confidence interval

		if self.globalMaximumExists() : 
			MLE = fminbound(self.help_loglikelihood, 0.0, 1.0, [], xtol=0.0000001) # TODO check! 
			max_logl = self.loglikelihood(MLE)
			
			self.MLE_cont = MLE 
			self.max_logl = max_logl 

			if self.loglikelihood(0.0) - max_logl + 1.92 >= 0: # 0.0 should be in the 95% confidence interval
				CI[0] = 0.0 
			else: # determine the zero point numerically - point should lie on the interval [0, MLE]
				CI[0] = bisect(self.help_ci_loglikelihood, 0.0, MLE)

			if self.loglikelihood(1.0) - max_logl + 1.92 >= 0: # 1.0 should be in the 95% confidence interval
				CI[1] = 1.0 
			else: # determine the zero point numerically - point should lie on the interval [MLE, 1]
				CI[1] = bisect(self.help_ci_loglikelihood, MLE, 1.0)
			return MLE, True, CI 

		else: 
			return None, False, CI

	def posteriorDistribution(self, vcf_record, paired_end_alignments, overlapping_alignments):
		"""
			Returns the posterior distribution for the vaf's in the vaf_range.
			Assumes an uniform prior over the vaf range.
		"""
		delta = getDelta (vcf_record)
		self.set(paired_end_alignments, overlapping_alignments, delta)
		
		posterior_probabilities = []
		sum_likelihoods = 0.0 
		
		for vaf in self.vaf_range: 
			p = math.exp(self.loglikelihood(vaf))
			posterior_probabilities.append(p)
			sum_likelihoods = sum_likelihoods + p
			
		if sum_likelihoods != 0.0: 
			return [p/sum_likelihoods for p in posterior_probabilities]
		else: 
			return None
		


class VAFEstimator:
	"""
		Performs variant allele frequency estimation on every indel in a 
		given VCF file on the basis of the alignments in BAM file
	"""
	def __init__(self, vcf_filename, bam_filename, mu, sigma, vaf_range=[0.0, 0.5, 1.0], search_range=5000):
		self.vcf_reader = vcf.Reader(open(vcf_filename))
		self.bam_processor = BAMProcessor(bam_filename, search_range)
		self.likelihood_maximizer = VAFLikelihoodMaximizer (mu, sigma, vaf_range = vaf_range)
		
		self.likelihood_maximizer.print()
		
		# TODO remove later
		self.test_data_collector = TestDataCollector()	

	def estimate(self): 
		# walk through all vcf records and estimate their VAF 
		for vcf_record in self.vcf_reader:
			if len(vcf_record.ALT) == 1: # records with several alternatives are ignored 
				if isDeletion (vcf_record) and returnIndelLength(vcf_record) >= 10: 
					deletion = Deletion(vcf_record)
					# Obtain relevant alignment data from the BAM processor 
					paired_end_alignments, overlapping_alignments = self.bam_processor.processDeletion(deletion)
					MLE_cont, unique_cont, CI = self.likelihood_maximizer.computeMLE_95CI(vcf_record, paired_end_alignments, overlapping_alignments)					
					MLE_disc, unique_disc = self.likelihood_maximizer.maximize(vcf_record, paired_end_alignments, overlapping_alignments)
					posterior_distr = self.likelihood_maximizer.posteriorDistribution(vcf_record, paired_end_alignments, overlapping_alignments)
					self.test_data_collector.updateDeletion(MLE_disc, unique_disc, MLE_cont, unique_cont, CI, deletion, paired_end_alignments, overlapping_alignments, posterior_distr) # TODO remove later
				#if isInsertion (vcf_record): 
				#	insertion = Insertion(vcf_record)
					# Obtain relevant alignment data from the BAM processor 
				#	paired_end_alignments, overlapping_alignments = self.bam_processor.processInsertion(insertion)
				#	MLE, unique = self.likelihood_maximizer.maximize(vcf_record, paired_end_alignments, overlapping_alignments)
				#	posterior_distr = self.likelihood_maximizer.posteriorDistribution(vcf_record, paired_end_alignments, overlapping_alignments)
					#self.test_data_collector.updateInsertion(MLE, unique, insertion, paired_end_alignments, overlapping_alignments, posterior_distr) # TODO remove later
					

		# After all is done, close the BAM file
		self.bam_processor.close()
		
		print("*********************************** CONFIDENCE LEVEN >= 0.0 ***********************************") 
		self.test_data_collector.summarize(confidence_level=0.0) # TODO remove later
		print("\n\n*********************************** CONFIDENCE LEVEN >= 0.9 ***********************************") 
		self.test_data_collector.summarize(confidence_level=0.9) # TODO remove later
		print("\n\n*********************************** CONFIDENCE LEVEN >= 0.99 ***********************************") 
		self.test_data_collector.summarize(confidence_level=0.99) # TODO remove later
	
					

				
def main():

	parser = OptionParser(usage=usage)
	
	parser.add_option("--unit", action="store_true", dest="continuous", default=False,
						help="The variant allele frequency (VAF) can take any value in the unit interval. Overwrites the '-p' option. (Default = false)")
	parser.add_option("-p", action="store", dest="ploidy", default=2, type=int, 
						help="Ploidy: the number of sets of chromosomes in the nucleus of a cell. (Default = 2; diploid)")
	parser.add_option("-c", action="store", dest="chromosome", default="", type=str, 
					help="In case one wants to analyze specifically one chromosome, e.g., '-c = 22' would entail that only indels on chr. 22 are consdired. (Default: all chromosomes)")
	parser.add_option("-r", action="store", dest="search_range", default=5000, type=int, 
					help="Range to search for potentially relevant reads (Default: 5000 bp)")
	(options, args) = parser.parse_args()
	
	if (len(args)!=4):
		parser.print_help()
		return 1

	# mean and standard deviation of the normal distribution that the internal segment lengths are assumed to follow when
	# not affected by an indel
	mu 	= float(args[0])
	sigma 	= float(args[1])
	
	vcf_filename = os.path.abspath(args[2])
	bam_filename = os.path.abspath(args[3])
	
	# Set the support for the variant allele frequency
	vaf_range = [] # VAF can vary over the unit interval	
	if not options.continuous: 
		vaf_range = [x/float(options.ploidy) for x in range(0,options.ploidy+1)] 

	vaf_estimator = VAFEstimator(vcf_filename, bam_filename, mu, sigma, vaf_range=vaf_range, search_range=options.search_range)
	vaf_estimator.estimate()

if __name__ == '__main__':
	sys.exit(main())
