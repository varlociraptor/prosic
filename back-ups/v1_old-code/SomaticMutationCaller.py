#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
from collections import defaultdict
from scipy.stats import norm 
from scipy.optimize import fminbound
from scipy.optimize import newton
import os
import sys
import vcf
import pysam
import math

# import from VAFEstimator.py
from VAFEstimator import *


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
	<bam-healthy>	BAM file containing the 
				alignments of the disease sample. 
				NOTE: File needs to be sorted 
				and indexed. 
	<bam-cancer>	BAM file containing the 
				alignments of the disease sample. 
				NOTE: File needs to be sorted 
				and indexed. 


NOTE: the current implementation
		1) 	ignores indels in the vcf file 
			for which more than one 
			alternative (ALT) are provided. 	
		2)	only considers primary 
			alignments; secondary/ternary etc. 
			alignments are neglected. 
"""


def ln(x):
	"""Natural logarithm. Returns -inf when x is equal to 0"""
	if x < sys.float_info.min:
		return float("-inf")
	return math.log(x)


class SMLikelihoodMaximizer: 
	"""
		Class that contains all functionality related to the (log)likelihood function
		for healthy and disease/control read data 
	"""
	def __init__(self, mu, sigma, alpha = 0.0, cancer_vaf_range = None, no_uncertainty = False):
		self.mu 	= float(mu)
		self.sigma 	= float(sigma)
		self.alpha 	= alpha
		self.cancer_vaf_range = cancer_vaf_range
		# Data for healthy sample
		self.h_paired_end_alignments = []	
		self.h_overlapping_alignments = []
		# Data for disease sample (cancer)
		self.c_paired_end_alignments = []
		self.c_overlapping_alignments = []

		self.delta = 0.0 
		self.epsilon0 = 0.0 
		self.epsilon1 = 0.0

		self.no_uncertainty = no_uncertainty # When True, every alignment is taken 100% serious 

	def setCancerVAFRange (self, cancer_vaf_range):
		"""
			Sets the cancer VAF range. If cancer_vaf_range = [] or 'None', the cancer_vaf_range is considered lie within 
			the unit interval. The MLE can, thus, range continuosly over that region.
		"""
		self.cancer_vaf_range = cancer_vaf_range

	def f(self, x, mean, std):
		"""Returns probability of observing x given a discretized normal distribution with mean mu and std of sigma, truncated below 0"""
		return (norm.cdf((x + 1.0 - mean) / std) - norm.cdf((x - mean) / std)) / (1.0 - norm.cdf(-1 * mean / std))


	def set(self, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments, delta):
		self.h_paired_end_alignments = h_paired_end_alignments
		self.h_overlapping_alignments = h_overlapping_alignments
		self.c_paired_end_alignments = c_paired_end_alignments
		self.c_overlapping_alignments = c_overlapping_alignments
		self.delta = float(delta)
		self.epsilon0 = returnEpsilon0(delta)
		self.epsilon1 = returnEpsilon1(delta)

		if self.no_uncertainty: # every read is taken 100% serious...
			print("Uncertainty is off!")
			for a in h_paired_end_alignments:
				a.probability = 1.0 
			for a in h_overlapping_alignments:
				a.probability = 1.0 
			for a in c_paired_end_alignments:
				a.probability = 1.0 
			for a in c_overlapping_alignments:
				a.probability = 1.0 

			self.epsilon0 = 0.0
			self.epsilon1 = 0.0




	def maximize(self, vcf_record, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments):
		delta = getDelta (vcf_record)
		
		if self.cancer_vaf_range == None or self.cancer_vaf_range == []: # cancer vaf is continuous and lies in the unit interval
			return self.computeMLE_continuous (h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments, delta)
		else:
			return self.computeMLE_discrete (h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments, delta)

	def loglikelihood(self, h_vaf, c_vaf):
		"""
			The loglikelihood function. h_vaf is the VAF of the healthy cells. c_vaf of the cancer cells
		"""
		logl = 0.0
		# First, we process the alignments of the control sample
		for paired_end_alignment in self.h_paired_end_alignments: # internal segment based evidence
			x = paired_end_alignment.value
			p = paired_end_alignment.probability
			logl = logl + ln(p * (h_vaf * self.f(x, self.mu + self.delta, self.sigma) + (1 - h_vaf) * self.f(x, self.mu, self.sigma)) + (1.0 - p))
		for overlapping_alignment in self.h_overlapping_alignments: # overlapping read evidence
			y = overlapping_alignment.value
			p = overlapping_alignment.probability 
			if y == 1: 
				logl = logl + ln(p * (h_vaf * (1.0 - self.epsilon1) + (1.0 - h_vaf) * self.epsilon0) + (1.0 - p))
			else: 
				logl = logl + ln(p * (h_vaf * self.epsilon1 + (1.0 - h_vaf) * (1.0 - self.epsilon0)) + (1.0 - p))


		# Secondly, we process the alignments of the disease/cancer sample
		for paired_end_alignment in self.c_paired_end_alignments: # internal segment based evidence
			x = paired_end_alignment.value
			p = paired_end_alignment.probability
			f_delta = self.f(x, self.mu + self.delta, self.sigma)
			f_0 = self.f(x, self.mu, self.sigma)
			logl = logl + ln(p * (self.alpha * (h_vaf * f_delta + (1 - h_vaf) * f_0) + (1.0 - self.alpha) * (c_vaf * f_delta + (1 - c_vaf) * f_0)) + (1.0 - p))
		for overlapping_alignment in self.c_overlapping_alignments: # overlapping read evidence
			y = overlapping_alignment.value
			p = overlapping_alignment.probability 
			if y == 1: 
				logl = logl + ln(p * (self.alpha * (h_vaf * (1.0 - self.epsilon1) + (1.0 - h_vaf) * self.epsilon0) + (1.0 - self.alpha) * (c_vaf * (1.0 - self.epsilon1) + (1.0 - c_vaf) * self.epsilon0)) + (1.0 - p))
			else: 
				logl = logl + ln(p * (self.alpha * (h_vaf * self.epsilon1 + (1.0 - h_vaf) * (1.0 - self.epsilon0)) + (1.0 - self.alpha) * (c_vaf * self.epsilon1 + (1.0 - c_vaf) * (1.0 - self.epsilon0))) + (1.0 - p))

		return logl 
		

	def computeMLE_discrete (self, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments, delta):
		"""
			Computes the MLE when the support of the VAF in the disease/cancer is discrete. 
			Returns both the MLE and a boolean, which is true when the MLE is unique and false otherwise.
		"""	
		self.set(h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments, delta)

		MLE = [[0.0, self.cancer_vaf_range[0]]]
		max_loglikelihood = self.loglikelihood (0.0, self.cancer_vaf_range[0]) 

		for i in range(1, len(self.cancer_vaf_range)):
			c_vaf = self.cancer_vaf_range[i]
			logl = self.loglikelihood(0.0, c_vaf)
			if logl > max_loglikelihood:
				max_loglikelihood = logl 
				MLE = [[0.0, c_vaf]]
			elif logl == max_loglikelihood:
				MLE.append([0.0, c_vaf])	

		for h_vaf in [0.5, 1.0]:
			for c_vaf in self.cancer_vaf_range:
				logl = self.loglikelihood(h_vaf, c_vaf)
				if logl > max_loglikelihood:
					max_loglikelihood = logl 
					MLE = [[h_vaf, c_vaf]]
				elif logl == max_loglikelihood:
					MLE.append([h_vaf, c_vaf])	
	
		if len(MLE) == 1: 
			return MLE[0], True
		else: 
			return MLE, False
	





	def globalMaximumExists (self):
		"""Determines whether a global maximum exists (for the continuous case)"""
		if self.alpha == 1.0: # no cancer cells present 
			return False

		for a in self.c_paired_end_alignments: 
			if a.probability > 0.0 and (self.f(a.value, self.mu + self.delta, self.sigma) is not self.f(a.value, self.mu, self.sigma)): 
				return True
		for a in self.c_overlapping_alignments: 
			if a.probability > 0.0: 
				return True
		return False

	def help_loglikelihood0(self, cancer_vaf):
		return -1 * self.loglikelihood(0.0, cancer_vaf)

	def help_loglikelihood1(self, cancer_vaf):
		return -1 * self.loglikelihood(0.5, cancer_vaf)

	def help_loglikelihood2(self, cancer_vaf):
		return -1 * self.loglikelihood(1.0, cancer_vaf)
			
	def computeMLE_continuous (self, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments, delta):
		"""
			Computes the MLE when the support of the VAF in the disease/cancer sample is the unit interval. 
			Returns both the MLE and a boolean, which is true when the 
			MLE is unique and false otherwise.
		"""
		self.set(h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments, delta)
		if self.globalMaximumExists() : 
			MLE0 = fminbound(self.help_loglikelihood0, 0.0, 1.0, [], xtol=0.0000001) # MLE for h_vaf = 0.0
			MLE1 = fminbound(self.help_loglikelihood1, 0.0, 1.0, [], xtol=0.0000001) # MLE for h_vaf = 0.5 
			MLE2 = fminbound(self.help_loglikelihood2, 0.0, 1.0, [], xtol=0.0000001) # MLE for h_vaf = 1.0 
			fval0 = self.loglikelihood(0.0, MLE0) 
			fval1 = self.loglikelihood(0.5, MLE1) 
			fval2 = self.loglikelihood(1.0, MLE2) 
			if fval0 > fval1 and fval0 > fval2:
				return [0.0, MLE0], fval0, True 
			elif fval2 > fval0 and fval2 > fval1:
				return [1.0, MLE2], fval2, True
			else:
				return [0.5, MLE1], fval1, True
		else: 
			return [None, None], None, False


	def help_ci_loglikelihood(self, c_vaf):
		return self.loglikelihood(self.mle_h_vaf, c_vaf) - self.max_logl + 1.92


	def computeMLE_95CI (self, vcf_record, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments):
		"""
			Computes the MLE when the support of the VAF in the disease/cancer sample is assumed to be the unit interval.
			In addition, the approx. 95% confidence interval is given. 
			Output:
				- Maximum likelihood estimate
				- Boolean; true, when MLE is unique and false otherwise
				- 95% confidence interval
		"""
		delta = getDelta (vcf_record)
		[mle_h_vaf, mle_c_vaf], max_logl, unique = self.computeMLE_continuous (h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments, delta)

		self.mle_h_vaf = mle_h_vaf
		self.max_logl = max_logl

		CI = [0.0, 1.0]

		if unique: 
			if self.loglikelihood(mle_h_vaf, 0.0) - max_logl + 1.92 >= 0: # 0 is in the 95% CI
				CI[0] = 0.0 
			else: # need to determine the zero point numerically
				CI[0] = bisect(self.help_ci_loglikelihood, 0.0, mle_c_vaf)
			
			if self.loglikelihood(mle_h_vaf, 1.0) - max_logl + 1.92 >= 0: # 1 is in the 95% CI
				CI[1] = 1.0 
			else: # need to determine the zero point numerically
				CI[1] = bisect(self.help_ci_loglikelihood, mle_c_vaf, 1.0)
			
			return [mle_h_vaf, mle_c_vaf], True, CI
				
		return [None, None], False, CI


	def probabilitySomaticMutation(self, vcf_record, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments):
		"""
			Returns the posterior probability of a somatic mutation. Only applies when the support of the VAF of the cancer cells is discrete.
		"""
		delta = getDelta(vcf_record)
		self.set(h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments, delta)

		prob_somatic_mutation = 0.0 
		prob_no_somatic_mutation = 0.0 

		for i in range(1, len(self.cancer_vaf_range)):
			c_vaf = cancer_vaf_range[i]
			prob_somatic_mutation = prob_somatic_mutation + math.exp(self.loglikelihood(0.0, c_vaf))

		prob_no_somatic_mutation = prob_no_somatic_mutation + math.exp(self.loglikelihood(0.0, 0.0))
		
		for h_vaf in [0.5, 1.0]:
			for c_vaf in self.cancer_vaf_range:
				prob_no_somatic_mutation = prob_no_somatic_mutation + math.exp(self.loglikelihood(h_vaf, c_vaf))
		
		return 	prob_somatic_mutation / (prob_somatic_mutation + prob_no_somatic_mutation)
			

	def posteriorDistribution(self, vcf_record, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments): 
		"""
			Returns the posterior distribution when the cancer's VAF has a discrete support
		"""
		delta = getDelta(vcf_record)
		self.set(h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments, delta)

		posterior_probabilities = [] 
		sum_likelihoods = 0.0 

		for h_vaf in [0.0, 0.5, 1.0]:
			post = []
			for c_vaf in self.cancer_vaf_range:
				p = math.exp(self.loglikelihood(h_vaf, c_vaf))
				post.append(p)
				sum_likelihoods = sum_likelihoods + p 
			posterior_probabilities.append(post)

		if sum_likelihoods != 0.0:
			for i in [0,1,2]:
				posterior_probabilities[i] = [p/sum_likelihoods for p in posterior_probabilities[i]]
			return posterior_probabilities 
		else:
			return None



class SMTestCase:
	"""
		Class that contains the relevant information for one indel.
		Used for testing the somatic mutation caller.

		TODO: remove for final version
	"""

	def __init__(self, is_deletion, h_true_vaf, c_true_vaf, true_sm, call, h_mle, c_mle, ci, unique, indel, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments):
		self.is_deletion 		= is_deletion
		self.h_true_vaf 		= h_true_vaf
		self.c_true_vaf			= c_true_vaf
		self.true_sm			= true_sm
		self.call 			= call 		
		self.h_mle 			= h_mle
		self.c_mle			= c_mle
		self.ci				= ci 
		self.unique			= unique
		self.indel			= indel
		self.position			= indel.vcf_record.POS
		self.length			= indel.length
		self.n_h_paired_end_alignments 	= len(h_paired_end_alignments)
		self.n_h_overlapping_alignments = len(h_overlapping_alignments)
		self.n_c_paired_end_alignments 	= len(c_paired_end_alignments)
		self.n_c_overlapping_alignments = len(c_overlapping_alignments)

		if self.true_sm == self.call:
			self.correct = True
		else:
			self.correct = False 
		
		self.Xh, self.Yh, self.Xc, self.Yc = [], [], [], []
	
		for a in h_paired_end_alignments: 
			self.Xh.append(a.value)

		for a in h_overlapping_alignments:
			self.Yh.append(a.value)

		for a in c_paired_end_alignments: 
			self.Xc.append(a.value)

		for a in c_overlapping_alignments:
			self.Yc.append(a.value)

	def print(self, min_length = 0):
		"""Prints relevant information to the command line"""
		if self.length >= min_length: 
			print(self.length, '\t|\t', self.correct, '\t|\t', self.true_sm, '\t',self.call, '\t|\t', self.h_true_vaf, self.c_true_vaf, '\t|\t', self.h_mle, self.c_mle, self.ci)

	def isOfLength(self, a, b):
		"""Determines whether this test case falls into a certain length class"""
		if self.length >= a and self.length <= b:
			return True
		return False


class SMTestDataCollector: 
	"""
		Collects/stores data used for testing

		TODO: remove for online/final version
	"""
	def __init__(self):
		self.testcases = [] 
		self.length_range = [[10, 29],[30,49],[50,99],[100,249]]
		self.min_length = 10 


	def determineVAF (self, gt_nums):
		if gt_nums == '1|1':
			return 1.0 
		elif gt_nums == None:
			return 0.0 
		return 0.5

	def getTruth (self, vcf_record):
		"""	Obtains the ground truth from the vcf_record: 
			- the true VAF for the healthy cells
			- the true cancer VAF  
			- somatic mutation or not
		"""
		h_vaf 		= 0.0  # healthy VAF
		# VAFs of the 4 cancer populations 
		som1_vaf 	= 0.0 
		som2_vaf 	= 0.0 
		som3_vaf 	= 0.0 
		som4_vaf 	= 0.0 

		for call in vcf_record.samples:
			if call.sample=='Control':
				h_vaf = self.determineVAF (call.gt_nums)
			if call.sample=='Som1':
				som1_vaf = self.determineVAF (call.gt_nums)	
			if call.sample=='Som2':
				som2_vaf = self.determineVAF (call.gt_nums)	
			if call.sample=='Som3':
				som3_vaf = self.determineVAF (call.gt_nums)
			if call.sample=='Som4':
				som4_vaf = self.determineVAF (call.gt_nums)
		
		cancer_vaf = (1.0 / 3.0) * som1_vaf + (1.0 / 3.0) * som2_vaf + (1.0 / 6.0) * som3_vaf + (1.0 / 6.0) * som4_vaf 

		somatic_mutation = False
		if h_vaf == 0.0 and cancer_vaf > 0.0:
			somatic_mutation = True

		return h_vaf, cancer_vaf, somatic_mutation  



	def updateDeletion (self, vcf_record, somatic_mutation, h_mle, c_mle, ci, unique, indel, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments):
		h_true_vaf, c_true_vaf, true_sm = self.getTruth(vcf_record)
		test_case = SMTestCase(True, h_true_vaf, c_true_vaf, true_sm, somatic_mutation, h_mle, c_mle, ci, unique, indel, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments)
		self.testcases.append(test_case)
		self.testcases[-1].print()

	def updateInsertion (self, vcf_record, somatic_mutation, h_mle, c_mle, ci, unique, indel, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments):
		h_true_vaf, c_true_vaf, true_sm = self.getTruth(vcf_record)
		test_case = SMTestCase(False, h_true_vaf, c_true_vaf, true_sm, somatic_mutation, h_mle, c_mle, ci, unique, indel, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments)
		self.testcases.append(test_case)
	
	def printTable2x2 (self, table, title):
		row_total, column_total = [0,0], [0,0]
		n = 0

		for r in [0,1]:
			for c in [0,1]:
				row_total[r] = row_total[r] + table[r][c]
				column_total[c] = column_total[c] + table[r][c]
				n = n + table[r][c]

		recall = None
		precision = None
		if column_total[0] != 0:
			recall = table[0][0] / float(column_total[0])
		if row_total[0] != 0:
			precision = table[0][0] / float(row_total[0])

		print('------------------------TABLE', title, '------------------------')
		print('\t\t\t\t\ttruth')
		print('\t\t|\tsomatic\t\t|\tnot somatic\t|\ttotal')
		print('--------------------------------------------------------------------------------------------')
		print('somatic\t\t|\t', table[0][0], '\t\t|\t', table[0][1], '\t\t|\t', row_total[0])
		print('not somatic\t|\t', table[1][0], '\t\t|\t', table[1][1], '\t\t|\t', row_total[1])	
		print('--------------------------------------------------------------------------------------------')
		print('total\t\t|\t', column_total[0], '\t\t|\t', column_total[1], '\t\t|\t', n)

		print('\nRECALL:\t\t', recall)
		print('PRECISION:\t', precision)	
		print('-------------------------------------------------------------------\n\n')	


	def summarize (self):
		"""Processes the gathered information, generates a summary and outputs it to the command line"""
		for range in self.length_range:
			total = 0.0			# total number of indels in this range 
			n_unique = 0.0			# number of indels for which the MLE was unique -> there was data to say something about their VAF
			
			table = [[0,0],[0,0]]		
			for test_case in self.testcases:
				if test_case.isOfLength(range[0], range[1]):
					total = total + 1
					if test_case.unique: 
						n_unique = n_unique + 1
						if test_case.true_sm and test_case.call:
							table[0][0] = table[0][0] + 1
						elif test_case.true_sm and not test_case.call:
							table[1][0] = table[1][0] + 1
						elif not test_case.true_sm and test_case.call:
							table[0][1] = table[0][1] + 1
						else: 
							table[1][1] = table[1][1] + 1

			print('*************************** LENGTH RANGE', range, '*******************************')
			print('total:\t\t', int(total))
			print('# unique:\t', int(n_unique), '\n')
			self.printTable2x2(table, 'Performance')
			print('*******************************************************************************************************')
	

	def printBoolean (self, boolean, f):
		if boolean:
			print('1\t', end = "", file = f)
		else:
			print('0\t', end = "", file = f)

	def printNumber (self, x, f):
		print(str(x) + '\t', end = "", file = f)

	def printList (self, l, f):
		s = '['
		if len(l) == 0:
			s = '[]'	
		else:
			for i in range(len(l) - 1):
				s = s + str(l[i]) + ', '
			s = s + str(l[-1]) + ']'
		print(s + '\t', end = "", file = f)

	def printToFile (self, output_filename):
		"""
			Stores the collected results into a file. 
		"""		
		f = open(output_filename, 'w')

		print('del\tchr\tpos\tlength\tcorrect\tsomatic\th_vaf\tc_vaf\tunique\tcall\th_mle\tc_mle\tci_start\tci_end\tX_h\tY_h\tX_c\tY_c', file = f) # header

		for t in self.testcases: 
			self.printBoolean(t.is_deletion, f)
			print(str(t.indel.chromosome) + '\t', end = "", file = f)
			self.printNumber(t.position, f)
			self.printNumber(t.length, f)
			self.printBoolean(t.correct, f)
			self.printBoolean(t.true_sm, f)
			self.printNumber(t.h_true_vaf, f)
			self.printNumber(t.c_true_vaf, f)
			self.printBoolean(t.unique, f)
			self.printBoolean(t.call, f)
			self.printNumber(t.h_mle, f)
			self.printNumber(t.c_mle, f)
			self.printNumber(t.ci[0], f)
			self.printNumber(t.ci[1], f)
			self.printList(t.Xh, f)
			self.printList(t.Yh, f)
			self.printList(t.Xc, f)	
			self.printList(t.Yc, f)
			print('\n', end="", file = f) 
		f.close()

class SomaticMutationCaller:
	
	def __init__(self, vcf_filename, bam_healthy_filename, bam_cancer_filename, mu, sigma, alpha = 0.0, cancer_vaf_range = None, search_range=5000, output_filename = "somatic_mutation_caller_output.txt", no_uncertainty = False):
		self.vcf_reader = vcf.Reader(open(vcf_filename))	
		self.bam_healthy_processor = BAMProcessor(bam_healthy_filename, search_range)	# BAM file for healthy (diploid) sample	
		self.bam_cancer_processor = BAMProcessor(bam_cancer_filename, search_range)	# BAM file for cancer/disease sample
		self.sm_likelihood_maximizer = SMLikelihoodMaximizer (mu, sigma, alpha = alpha, cancer_vaf_range = cancer_vaf_range, no_uncertainty = no_uncertainty)		# Maximizes the likelihood 

		# TODO remove later on
		self.sm_test_data_collector = SMTestDataCollector()
		self.output_filename = output_filename

	def somaticMutation(self, MLE_healthy_VAF, CI_cancer_VAF):
		if MLE_healthy_VAF == 0.0 and CI_cancer_VAF[0] > 0.0:
			return True
		else:
			return False

	def call(self):
		"""
			Classifies the indels present in the vcf file as somatic/not somatic given the BAM alignments in the healthy and cancer/disease sample 
		"""
		#i = 0 

		# walk through all vcf records that represent indels and call them as being somatic or not
		for vcf_record in self.vcf_reader: 
			if len(vcf_record.ALT) == 1 and vcf_record.CHROM == '1': # records with several alternatives are ignored 
				if isDeletion (vcf_record) and returnIndelLength(vcf_record) >= 10: 
					#i = i + 1
					print('Deletion @', vcf_record.CHROM, ' pos.', vcf_record.POS)
					deletion = Deletion(vcf_record)
					
					# Obtain relevant alignment data from the control (h) and disease/cancer (c) BAM files
					h_paired_end_alignments, h_overlapping_alignments = self.bam_healthy_processor.processDeletion (deletion) # healthy 
					c_paired_end_alignments, c_overlapping_alignments = self.bam_cancer_processor.processDeletion (deletion) # cancer

					[MLE_healthy_VAF, MLE_cancer_VAF], unique, CI_cancer_VAF = self.sm_likelihood_maximizer.computeMLE_95CI (vcf_record, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments)
					# Perform the actual call on the basis of the MLEs and the 95% confidence interval (CI)
					somatic_mutation = self.somaticMutation(MLE_healthy_VAF, CI_cancer_VAF)

					self.sm_test_data_collector.updateDeletion(vcf_record, somatic_mutation, MLE_healthy_VAF, MLE_cancer_VAF, CI_cancer_VAF, unique, deletion, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments) 
					

					#posterior_distr = self.sm_likelihood_maximizer.posteriorDistribution(vcf_record, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments)
					#self.test_data_collector.updateDeletion(MLE, unique, deletion, paired_end_alignments, overlapping_alignments, posterior_distr) # TODO remove later
				#if isInsertion (vcf_record): 
				#	insertion = Insertion(vcf_record)
					# Obtain relevant alignment data from the disease and control BAM files  
				#	h_paired_end_alignments, h_overlapping_alignments = self.bam_healthy_processor.processInsertion (insertion) # healthy 
				#	c_paired_end_alignments, c_overlapping_alignments = self.bam_cancer_processor.processInsertion (insertion) # cancer

					#somatic_mutation, MLE, unique = self.sm_likelihood_maximizer.maximize(vcf_record, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments)

					#posterior_distr = self.sm_likelihood_maximizer.posteriorDistribution(vcf_record, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments)
					#self.test_data_collector.updateInsertion(MLE, unique, insertion, paired_end_alignments, overlapping_alignments, posterior_distr) # TODO remove later
			
			#if i == 1000:
				#break  

		self.sm_test_data_collector.summarize()		
		self.sm_test_data_collector.printToFile(self.output_filename)
		
		# After all is done, close the BAM files
		self.bam_healthy_processor.close()
		self.bam_cancer_processor.close()




def main():


	parser = OptionParser(usage=usage)
	
	parser.add_option("-a", action="store", dest="alpha", default=0.0, type=float, 
					help="Level of impurity in disease/cancer sample. (Default = 0.0; no healthy cells present in the disease/cancer sample)")
	parser.add_option("-p", action="store", dest="ploidy", default=None, type=int, 
						help="Ploidy: the number of sets of chromosomes in the nucleus of a cancer cell. (Default = no assumptions made; unit interval is used)")
	parser.add_option("-c", action="store", dest="chromosome", default="", type=str, 
					help="In case one wants to analyze specifically one chromosome, e.g., '-c = 22' would entail that only indels on chr. 22 are consdired. (Default: all chromosomes)")
	parser.add_option("-r", action="store", dest="search_range", default=5000, type=int, 
					help="Range to search for potentially relevant reads (Default: 5000 bp)")
	parser.add_option("--uncertainty_off", action="store_true", dest="uncertainty_off", default=False,
						help="Alignment uncertainties are neglected. (Default = false)")
	parser.add_option("-o", action="store", dest="output_file", default="somatic_mutation_caller_output.txt", 
						help="Output file for the results. Only used for testing. (Default = somatic_mutation_caller_output.txt)")
	(options, args) = parser.parse_args()

	if (len(args)!=5):
		parser.print_help()
		return 1

	# mean and standard deviation of the normal distribution that the internal segment lengths are assumed to follow when
	# not affected by an indel
	mu 	= float(args[0])
	sigma 	= float(args[1])
	
	vcf_filename = os.path.abspath(args[2])
	bam_healthy_filename = os.path.abspath(args[3])
	bam_cancer_filename = os.path.abspath(args[4])
	
	# Set the support for the variant allele frequency
	cancer_vaf_range = [] # VAF can vary over the unit interval	
	if options.ploidy is not None: 
		cancer_vaf_range = [x/float(options.ploidy) for x in range(0,options.ploidy+1)] 
	else:
		cancer_vaf_range = None

	sm_caller = SomaticMutationCaller(vcf_filename, bam_healthy_filename, bam_cancer_filename, mu, sigma, alpha = options.alpha, cancer_vaf_range = cancer_vaf_range, search_range = options.search_range, output_filename = options.output_file, no_uncertainty = options.uncertainty_off) 
	sm_caller.call() 

if __name__ == '__main__':
	sys.exit(main())
