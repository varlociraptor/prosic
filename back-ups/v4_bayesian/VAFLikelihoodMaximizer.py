#!/usr/bin/env python
from __future__ import print_function, division
from scipy.stats import norm 
from scipy.optimize import fminbound
from scipy.optimize import bisect
import sys
import math

from Indel import *
from Alignments import *

__author__ = "Louis Dijkstra"

"""
	VAFLikelihoodMaxiizer.py contains the functionality related to the VAF likelihood function. 
	Used mainly by the VAF Estimator 
"""

def ln(x):
	"""Natural logarithm. Returns -inf when x is equal to 0 (or below system precision)"""
	if x < sys.float_info.min:
		return float("-inf")
	return math.log(x)

class VAFLikelihoodMaximizer: 
	"""Class containing the functionality related to the VAF likelihood function."""

	def __init__(self, mu, sigma, epsilon0 = 0.0, epsilon1 = 0.0, vaf_range = [0.0, 0.5, 1.0], tolerance=0.0000001):
		self.mu 			= float(mu) 
		self.sigma 			= float(sigma)
		self.epsilon0 			= epsilon0 
		self.epsilon1 			= epsilon1
		self.vaf_range 			= vaf_range # if [] or 'None', the VAF's support is considered to be continuous on the unit interval [0,1]
		self.tolerance			= tolerance # precision required for continuos MLE
		# initialization
		self.paired_end_alignments 	= []
		self.overlapping_alignments 	= []
		self.delta 			= 0.0
		
	def set(self, paired_end_alignments, overlapping_alignments, delta):
		"""Provides the likelihood maximizer with the available data."""
		self.paired_end_alignments 	= paired_end_alignments
		self.overlapping_alignments 	= overlapping_alignments
		self.delta 			= float(delta)


	def f(self, x, mean, std):
		"""Returns the probability of observing x given a discretized normal distribution with mean mu and std of sigma, truncated below 0"""
		return (norm.cdf((x + 1.0 - mean) / std) - norm.cdf((x - mean) / std)) / (1.0 - norm.cdf(-1 * mean / std))
	
	def loglikelihood (self, vaf):
		"""Returns the loglikelihood."""
		logl = 0.0 
		for a in self.paired_end_alignments: # internal segment based evidence
			logl = logl + ln(a.probability * (vaf * self.f(a.value, self.mu + self.delta, self.sigma) + (1 - vaf) * self.f(a.value, self.mu, self.sigma)) + (1.0 - a.probability))
		for a in self.overlapping_alignments: # overlapping read evidence
			p = a.probability 
			if a.value == 1: 
				logl = logl + ln(a.probability * (vaf * (1.0 - self.epsilon1) + (1.0 - vaf) * self.epsilon0) + (1.0 - a.probability)) 
			else: 
				logl = logl + ln(a.probability * (vaf * self.epsilon1 + (1.0 - vaf) * (1.0 - self.epsilon0)) + (1.0 - a.probability))
		return logl 

	def maximize (self, vcf_record, paired_end_alignments, overlapping_alignments):
		"""Returns the maximum likelihood estimate of the VAF given an indel (vcf_record) and internal segment and overlapping alignments."""
		if self.vaf_range == None or self.vaf_range == []: # VAF lies within [0,1] and is continuous
			return self.computeMLE_continuous (vcf_record, paired_end_alignments, overlapping_alignments)  
		else:
			return self.computeMLE_discrete (vcf_record, paired_end_alignments, overlapping_alignments)
			
	def computeMLE_discrete(self, vcf_record, paired_end_alignments, overlapping_alignments): 
		"""Computes the MLE when the support of the VAF is discrete. Returns both the MLE and a boolean, which is true when the 
		   MLE is unique and false otherwise."""
		self.set(paired_end_alignments, overlapping_alignments, getDelta(vcf_record))
		
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

	def posteriorDistribution(self, vcf_record, paired_end_alignments, overlapping_alignments):
		"""Returns the posterior distribution for the vaf's in the vaf_range. Assumes an uniform prior over the vaf range."""
		self.set(paired_end_alignments, overlapping_alignments, getDelta(vcf_record))
		
		posterior_probabilities = []
		sum_likelihoods = 0.0 
		
		for vaf in self.vaf_range: 
			p = math.exp(self.loglikelihood(vaf))
			posterior_probabilities.append(p)
			sum_likelihoods = sum_likelihoods + p
			
		if sum_likelihoods != 0.0: 
			return [p/sum_likelihoods for p in posterior_probabilities]
		return [1.0 / len(self.vaf_range)] * len(self.vaf_range) # no info -> return uniform prior distribution
	
	def globalMaximumExists (self):
		"""Determines whether a unique global maximum exists when the VAF is continuous."""
		for a in self.paired_end_alignments: 
			if a.probability > 0.0 and (self.f(a.value, self.mu + self.delta, self.sigma) is not self.f(a.value, self.mu, self.sigma)): 
				return True
		if self.epsilon0 == 0.5 and self.epsilon1 == 0.5:
			return False
		for a in self.overlapping_alignments: 
			if a.probability > 0.0: 
				return True
		return False
	
	def help_loglikelihood(self, vaf):
		"""Auxiliary function used for numerically maximizing the loglikelihood function."""
		return -1 * self.loglikelihood(vaf)
	
	def computeMLE_continuous (self, vcf_record, paired_end_alignments, overlapping_alignments):
		"""Computes the MLE when the support of the VAF is the unit interval. Returns both the MLE and a boolean, which is true when the 
		   MLE is unique and false otherwise."""
		self.set(paired_end_alignments, overlapping_alignments, getDelta(vcf_record))
		if self.globalMaximumExists() : 
			return fminbound(self.help_loglikelihood, 0.0, 1.0, [], xtol=self.tolerance), True 
		else: 
			return None, False
		
	def help_ci_loglikelihood(self, vaf):
		"""Auxiliary function for computing the 95% confidence interval."""
		return self.loglikelihood(vaf) - self.max_logl + 1.92

	def computeMLE_95CI (self, vcf_record, paired_end_alignments, overlapping_alignments):
		"""Computes the MLE when the support of the VAF is the unit interval. In addition, the approx. 95% confidence interval is given. 
		   Output:
				- Maximum likelihood estimate
				- Boolean; true, when MLE is unique and false otherwise
				- 95% confidence interval"""
		self.set(paired_end_alignments, overlapping_alignments, getDelta(vcf_record))

		if not self.globalMaximumExists():
			return None, False, [0.0,1.0]

		self.MLE_cont = fminbound(self.help_loglikelihood, 0.0, 1.0, [], xtol=self.tolerance) # Determine MLE
		self.max_logl = self.loglikelihood(self.MLE_cont)

		CI = [0.0, 1.0] # confidence interval

		if self.loglikelihood(0.0) - self.max_logl + 1.92 < 0: 
			CI[0] = bisect(self.help_ci_loglikelihood, 0.0, self.MLE_cont)
		if self.loglikelihood(1.0) - self.max_logl + 1.92 < 0: 
			CI[1] = bisect(self.help_ci_loglikelihood, self.MLE_cont, 1.0)
		return self.MLE_cont, True, CI 


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

