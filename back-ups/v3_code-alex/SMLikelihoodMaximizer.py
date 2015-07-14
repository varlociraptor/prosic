#!/usr/bin/env python
from __future__ import print_function, division
from scipy.stats import norm 
from scipy.stats import chi2 
from scipy.optimize import fminbound
from scipy.optimize import bisect

import sys
import math

from Indel import *
from Alignments import *

__author__ = "Louis Dijkstra"

"""
	SMLikelihoodMaxiizer.py contains the functionality related to the likelihood function for the somatic mutation calling setting. 
	Used mainly by the somatic mutation caller. 
"""

def ln(x):
	"""Natural logarithm. Returns -inf when x is equal to 0 (or below system precision)"""
	if x < sys.float_info.min:
		return float("-inf")
	return math.log(x)
	
class SMLikelihoodMaximizer: 
	"""Class containing all functionality related to the (log)likelihood function for somatic mutation calling""" 
	def __init__(self, mu, sigma, alpha = 0.0, epsilon0 = 0.0, epsilon1 = 0.0, cancer_vaf_range = None, tolerance=0.0000001):
		self.mu 		= float(mu)		
		self.sigma 		= float(sigma)
		self.alpha 		= alpha
		self.epsilon0		= epsilon0 
		self.epsilon1 		= epsilon1
		self.cancer_vaf_range 	= cancer_vaf_range # if [] or 'None', the cancer VAF's support is considered to be continuous on the unit interval [0,1]
		self.tolerance		= tolerance # precision required for continuos MLE

		# Initialization 
		# Data for healthy sample 
		self.h_paired_end_alignments 	= []	
		self.h_overlapping_alignments 	= []
		# Data for disease sample (cancer)
		self.c_paired_end_alignments 	= []
		self.c_overlapping_alignments 	= []

		self.delta 	= 0.0 
		self.epsilon0 	= 0.0 
		self.epsilon1 	= 0.0

	def set(self, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments, delta):
		"""Provides the likelihood maximizer with the available data."""
		self.h_paired_end_alignments 	= h_paired_end_alignments
		self.h_overlapping_alignments 	= h_overlapping_alignments
		self.c_paired_end_alignments 	= c_paired_end_alignments
		self.c_overlapping_alignments 	= c_overlapping_alignments
		self.delta 			= float(delta)

	def f(self, x, mean, std):
		"""Returns probability of observing x given a discretized normal distribution with mean mu and std of sigma, truncated below 0"""
		return (norm.cdf((x + 1.0 - mean) / std) - norm.cdf((x - mean) / std)) / (1.0 - norm.cdf(-1 * mean / std))

	def loglikelihood(self, h_vaf, c_vaf):
		"""Returns the loglikelihood. h_vaf and c_vaf denote the VAF of the healthy (h) and cancer (c) cells, respectively."""
		logl = 0.0
		# First, we process the alignments of the control sample
		for a in self.h_paired_end_alignments: # internal segment based evidence
			logl = logl + ln(a.probability * (h_vaf * self.f(a.value, self.mu + self.delta, self.sigma) + (1 - h_vaf) * self.f(a.value, self.mu, self.sigma)) + (1.0 - a.probability)) 
		for a in self.h_overlapping_alignments: # overlapping read evidence
			if a.value == 1: 
				logl = logl + ln(a.probability * (h_vaf * (1.0 - self.epsilon1) + (1.0 - h_vaf) * self.epsilon0) + (1.0 - a.probability)) 
			else: 
				logl = logl + ln(a.probability * (h_vaf * self.epsilon1 + (1.0 - h_vaf) * (1.0 - self.epsilon0)) + (1.0 - a.probability))

		# Secondly, we process the alignments of the disease/cancer sample
		for a in self.c_paired_end_alignments: # internal segment based evidence
			f_delta = self.f(a.value, self.mu + self.delta, self.sigma)
			f_0 	= self.f(a.value, self.mu, self.sigma)
			logl 	= logl + ln(a.probability * (self.alpha * (h_vaf * f_delta + (1.0 - h_vaf) * f_0) + (1.0 - self.alpha) * (c_vaf * f_delta + (1.0 - c_vaf) * f_0)) + (1.0 - a.probability)) 
		for a in self.c_overlapping_alignments: # overlapping read evidence
			if a.value == 1: 
				logl = logl + ln(a.probability * (self.alpha * (h_vaf * (1.0 - self.epsilon1) + (1.0 - h_vaf) * self.epsilon0) + (1.0 - self.alpha) * (c_vaf * (1.0 - self.epsilon1) + (1.0 - c_vaf) * self.epsilon0)) + (1.0 - a.probability)) 
			else: 
				logl = logl + ln(a.probability * (self.alpha * (h_vaf * self.epsilon1 + (1.0 - h_vaf) * (1.0 - self.epsilon0)) + (1.0 - self.alpha) * (c_vaf * self.epsilon1 + (1.0 - c_vaf) * (1.0 - self.epsilon0))) + (1.0 - a.probability)) 

		return logl 

	def maximize(self, vcf_record, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments):
		"""Returns the maximum likelihood estimate of the VAF of the healthy and cancer cells given an indel (vcf_record) and internal segment and overlapping alignments."""
		if self.cancer_vaf_range == None or self.cancer_vaf_range == []: # cancer vaf is continuous and lies in the unit interval
			return self.computeMLE_continuous (vcf_record, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments)
		else:
			return self.computeMLE_discrete (vcf_record, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments)

	def computeMLE_discrete (self, vcf_record, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments):
		"""Computes the MLE when the support of the cancer VAF is discrete. Returns both the MLE and a boolean, which is true when the 
		   MLE is unique and false otherwise."""
		self.set(h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments, getDelta(vcf_record))

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
	
		if len(MLE) == 1: # MLE is unique
			return MLE[0], True
		return MLE, False

	def posteriorDistribution(self, vcf_record, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments): 
		"""Returns the posterior distribution for the healthy VAF and the cancer VAF (both discrete). Assumes an uniform prior."""
		self.set(h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments, getDelta(vcf_record))

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
		else: # if the sum of the likelihoods is zero, we return the uniform prior distribution
			return [[1.0/(3 * len(self.cancer_vaf_range))] * len(self.cancer_vaf_range)] * 3

	def probabilitySomaticMutation(self, vcf_record, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments):
		"""Returns the posterior probability of a somatic mutation. Only applies when the support of the VAF of the cancer cells is discrete."""
		self.set(h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments, getDelta(vcf_record))

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
			
	def globalMaximumExists (self):
		"""Determines whether a unique global maximum exists for the cancer VAF when it is assumed continuous"""
		if self.alpha == 1.0: # no cancer cells present 
			return False
		for a in self.c_paired_end_alignments: # check whether there is information in the cancer cells 
			if a.probability > 0.0 and (self.f(a.value, self.mu + self.delta, self.sigma) is not self.f(a.value, self.mu, self.sigma)): 
				return True
		if self.epsilon0 == 0.5 and self.epsilon1 == 0.5:
			return False 
		for a in self.c_overlapping_alignments: 
			if a.probability > 0.0: 
				return True
		return False

	def help_loglikelihood0(self, cancer_vaf):
		"""Auxiliary function for numerically maximizing the cancer VAF when the healthy VAF is equal to 0.0"""
		return -1 * self.loglikelihood(0.0, cancer_vaf)

	def help_loglikelihood1(self, cancer_vaf):
		"""Auxiliary function for numerically maximizing the cancer VAF when the healthy VAF is equal to 0.5"""
		return -1 * self.loglikelihood(0.5, cancer_vaf)

	def help_loglikelihood2(self, cancer_vaf):
		"""Auxiliary function for numerically maximizing the cancer VAF when the healthy VAF is equal to 1.0"""
		return -1 * self.loglikelihood(1.0, cancer_vaf)
			
	def computeMLE_continuous (self, vcf_record, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments):
		"""Computes the MLE when the support of the cancer VAF is the unit interval (continuous). Returns both the MLE, the maximum loglikelihood and a boolean, which is true when the 
		   MLE is unique and false otherwise."""
		self.set(h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments, getDelta(vcf_record))

		if self.globalMaximumExists() : 
			MLE0 = fminbound(self.help_loglikelihood0, 0.0, 1.0, [], xtol=self.tolerance) # MLE for h_vaf = 0.0
			MLE1 = fminbound(self.help_loglikelihood1, 0.0, 1.0, [], xtol=self.tolerance) # MLE for h_vaf = 0.5 
			MLE2 = fminbound(self.help_loglikelihood2, 0.0, 1.0, [], xtol=self.tolerance) # MLE for h_vaf = 1.0 
			logl0 = self.loglikelihood(0.0, MLE0) 
			logl1 = self.loglikelihood(0.5, MLE1) 
			logl2 = self.loglikelihood(1.0, MLE2) 

			if logl0 == max([logl0, logl1, logl2]):
				return [0.0, MLE0], logl0, True	
			elif logl1 == max([logl0, logl1, logl2]):
				return [0.5, MLE1], logl1, True	
			else:
				return [1.0, MLE2], logl2, True	
		else: 
			return [None, None], None, False


	def returnPValue(self, c_vaf_mle):
		"""Returns the likelihood ratio test p-value. Only applies when the healthy VAF is 0."""
		T = 2.0*(self.loglikelihood(0.0, c_vaf_mle) - self.loglikelihood(0.0, 0.0))
		return 1.0 - chi2.cdf(T,1)

	def help_ci_loglikelihood(self, c_vaf):
		"""Auxiliary function for computing the 95% confidence interval of the cancer VAF."""
		return self.loglikelihood(self.mle_h_vaf, c_vaf) - self.max_logl + self.threshold

	def computeCI (self, h_vaf, mle_c_vaf, max_logl, confidence_level = 0.95):
		"""Returns the (confidence_level)*100 % confidence interval for the cancer VAF when the healthy VAF is equal to h_vaf."""
		CI = [0.0, 1.0]

		self.mle_h_vaf 	= h_vaf 
		self.max_logl 	= max_logl
		self.threshold 	= chi2.ppf(confidence_level, 1) / 2.0
		
		if self.loglikelihood(h_vaf, 0.0) - max_logl + self.threshold < 0:
			CI[0] = bisect(self.help_ci_loglikelihood, 0.0, mle_c_vaf)
		if self.loglikelihood(h_vaf, 1.0) - max_logl + self.threshold < 0:
			CI[1] = bisect(self.help_ci_loglikelihood, mle_c_vaf, 1.0) 
	
		return CI 

	def computeMLE_95CI (self, vcf_record, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments, confidence_level = 0.95):
		"""Computes the MLE when the support of the cancer VAF is the unit interval. In addition, the approx. (confidence_level)*100 % confidence interval is given. 
		   Output:
				- Maximum likelihood estimate
				- Boolean; true, when MLE is unique and false otherwise
				- (confidence_level)*100 % confidence interval"""
		[mle_h_vaf, mle_c_vaf], max_logl, unique = self.computeMLE_continuous (vcf_record, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments)
		
		if not unique:
			return [None, None], False, [0.0, 1.0]

		CI = self.computeCI(mle_h_vaf, mle_c_vaf, max_logl, confidence_level = confidence_level) 
			
		return [mle_h_vaf, mle_c_vaf], True, CI




