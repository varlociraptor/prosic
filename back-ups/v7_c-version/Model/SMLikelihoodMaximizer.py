#!/usr/bin/env python
from __future__ import print_function, division
from scipy.stats import norm 
from scipy.stats import chi2 
from scipy.optimize import fminbound
from scipy.optimize import bisect
from scipy.integrate import quad

import sys
import math
import numpy as np

from Variant.Indel import *
from BAM.Alignments import *

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

ln_np = np.vectorize(ln)

def isize_distr(isize, mean, std):
	"""Returns the probability of observing an insert size given the mean and standard deviation."""
	return (norm.cdf((isize + 1.0 - mean) / std) - norm.cdf((isize - mean) / std)) / (1.0 - norm.cdf(-1 * mean / std))

def isize_healthy_likelihood(isize, align_prob, h_vaf, mean, std, delta):
	"""Returns the likelihood of an individual insert size observation of the healthy sample."""	
	return align_prob * (h_vaf * isize_distr(isize, mean + delta, std) + (1.0 - h_vaf) * isize_distr(isize, mean, std)) + (1.0 - align_prob)

isize_healthy_likelihood_np = np.vectorize(isize_healthy_likelihood)	

def split_healthy_likelihood(split, align_prob, h_vaf, epsilon_a, epsilon_p):
	"""Returns the likelihood of an individual split read observation in the healthy sample."""
	return align_prob * (split * (h_vaf * (1.0 - epsilon_p) + (1.0 - h_vaf) * epsilon_a) + (1.0 - split) * (h_vaf * epsilon_p + (1.0 - h_vaf) * (1.0 - epsilon_a))) + (1.0 - align_prob) 

split_healthy_likelihood_np = np.vectorize(split_healthy_likelihood)	

def isize_cancer_likelihood(isize, align_prob, h_vaf, c_vaf, alpha, mean, std, delta):
	"""Returns the likelihood of an individual insert size observation of the healthy sample."""
	a = isize_distr(isize, mean, std)
	p = isize_distr(isize, mean + delta, std)	
	return align_prob * (alpha * (h_vaf * p + (1.0 - h_vaf) * a) + (1.0 - alpha)*(c_vaf * p + (1.0 - c_vaf)*a)) + + (1.0 - align_prob)

isize_cancer_likelihood_np = np.vectorize(isize_cancer_likelihood)

def split_cancer_likelihood(split, align_prob, h_vaf, c_vaf, alpha, epsilon_a, epsilon_p):
	"""Returns the likelihood of an individual split read observation in the healthy sample."""
	return align_prob * (alpha * (split * (h_vaf * (1.0 - epsilon_p) + (1.0 - h_vaf) * epsilon_a) + (1.0 - split) * (h_vaf * epsilon_p + (1.0 - h_vaf) * (1.0 - epsilon_a))) + (1.0 - alpha) * (split * (c_vaf * (1.0 - epsilon_p) + (1.0 - c_vaf) * epsilon_a) + (1.0 - split) * (c_vaf * epsilon_p + (1.0 - c_vaf) * (1.0 - epsilon_a)))) + (1.0 - align_prob) 

split_cancer_likelihood_np = np.vectorize(split_cancer_likelihood)	


class SMLikelihoodMaximizer:
	"""Class containing all functionality related to the (log)likelihood function for somatic mutation calling."""
	def __init__(self, mu_h, sigma_h, mu_c, sigma_c, alpha=0.0, epsilon_a=0.0, epsilon_p=0.0, cancer_vaf_range=None, tolerance=0.0000001):
		self.mu_h 		= float(mu_h)		
		self.sigma_h 		= float(sigma_h)
		self.mu_h 		= float(mu_c)		
		self.sigma_h 		= float(sigma_c)
		self.alpha 		= alpha
		self.epsilon_a		= epsilon_a 
		self.epsilon_p 		= epsilon_p
		self.cancer_vaf_range 	= cancer_vaf_range # if [] or 'None', the cancer VAF's support is considered to be continuous on the unit interval [0,1]
		self.tolerance		= tolerance # precision required for continuos MLE

	def likelihood(self, h_vaf, c_vaf, delta, h_isize, h_isize_prob, h_splits, h_splits_prob, c_isize, c_isize_prob, c_splits, c_splits_prob):
		"""Returns the likelihood. h_vaf and c_vaf denote the VAF of the healthy (h) and cancer (c) cells, respectively."""
		likelihood_isize_healthy 	= 1.0 if h_isize.size == 0 else isize_healthy_likelihood_np(h_isize, h_isize_prob, h_vaf, self.mu_h, self.sigma_h, delta)
		likelihood_splits_healthy 	= 1.0 if h_splits.size == 0 else split_healthy_likelihood_np(h_splits, h_splits_prob, h_vaf, self.epsilon_a, self.epsilon_p)
		likelihood_isize_cancer 	= 1.0 if c_isize.size == 0 else isize_cancer_likelihood_np(c_isize, c_isize_prob, h_vaf, c_vaf, self.alpha, self.mu_h, self.sigma_h, delta)
		likelihood_splits_cancer 	= 1.0 if c_splits.size == 0 else split_cancer_likelihood_np(c_splits, c_splits_prob, h_vaf, c_vaf, self.alpha, self.epsilon_a, self.epsilon_p)
		return np.prod(likelihood_isize_healthy) * np.prod(likelihood_splits_healthy) * np.prod(likelihood_isize_cancer) * np.prod(likelihood_splits_cancer)

	def loglikelihood(self, h_vaf, c_vaf, delta, h_isize, h_isize_prob, h_splits, h_splits_prob, c_isize, c_isize_prob, c_splits, c_splits_prob):
		logl_isize_healthy 	= 0.0 if h_isize.size == 0 else ln_np(isize_healthy_likelihood_np(h_isize, h_isize_prob, h_vaf, self.mu_h, self.sigma_h, delta))
		logl_splits_healthy 	= 0.0 if h_splits.size == 0 else ln_np(split_healthy_likelihood_np(h_splits, h_splits_prob, h_vaf, self.epsilon_a, self.epsilon_p))
		logl_isize_cancer 	= 0.0 if c_isize.size == 0 else ln_np(isize_cancer_likelihood_np(c_isize, c_isize_prob, h_vaf, c_vaf, self.alpha, self.mu_h, self.sigma_h, delta))
		logl_splits_cancer 	= 0.0 if c_splits.size == 0 else ln_np(split_cancer_likelihood_np(c_splits, c_splits_prob, h_vaf, c_vaf, self.alpha, self.epsilon_a, self.epsilon_p))
		return np.sum(logl_isize_healthy) + np.sum(logl_splits_healthy) + np.sum(logl_isize_cancer) + np.sum(logl_splits_cancer)

	def likelihood_aux(self, c_vaf, h_vaf, delta, h_isize, h_isize_prob, h_splits, h_splits_prob, c_isize, c_isize_prob, c_splits, c_splits_prob):
		return -1.0 * self.likelihood(h_vaf, c_vaf, delta, h_isize, h_isize_prob, h_splits, h_splits_prob, c_isize, c_isize_prob, c_splits, c_splits_prob)
	
	def loglikelihood_aux(self, c_vaf, h_vaf, delta, h_isize, h_isize_prob, h_splits, h_splits_prob, c_isize, c_isize_prob, c_splits, c_splits_prob):
		return -1.0 * self.loglikelihood(h_vaf, c_vaf, delta, h_isize, h_isize_prob, h_splits, h_splits_prob, c_isize, c_isize_prob, c_splits, c_splits_prob)

	def returnMLECancerVAF(self, h_vaf, delta, h_isize, h_isize_prob, h_splits, h_splits_prob, c_isize, c_isize_prob, c_splits, c_splits_prob):
		return fminbound(self.loglikelihood_aux, 0.0, 1.0, args=(h_vaf, delta, h_isize, h_isize_prob, h_splits, h_splits_prob, c_isize, c_isize_prob, c_splits, c_splits_prob), xtol=self.tolerance)

