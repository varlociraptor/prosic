#!/usr/bin/env python
from __future__ import print_function, division
from scipy.stats import norm 
from scipy.stats import chi2 

__author__ = "Louis Dijkstra"

"""
	ResultsInstance.py contains the functionailty to process intermediate results. 
"""

def obtainResults(results_file, all_chromosomes = False, confidence_level = 0.0, threshold_internal_segment_evidences = 0.0, threshold_split_read_evidences = 0.0):
	"""Returns all the results from a results file from chromosome 1. When the boolean all_chromosomes = True, then all results are returned."""
	results 	= [] # list with all the results
	true_c_vafs 	= [] # a list of all the true different cancer VAFs present in the file

	truth_vcf_file_used = False
	if 'truth' in results_file.name:
		truth_vcf_file_used = True
	
	# walk through all the results in the file
	for line in results_file: 

		chrom = line[3:5]
		if truth_vcf_file_used:
			chrom = line[:2]

		if chrom == 'X' or chrom == 'X ' or chrom == 'Y' or chrom == 'Y ' or chrom == '': # sex chromosomes are not taken into account
				continue
 
		if not all_chromosomes:			
			if int(chrom) > 3:
				continue 

		result = ResultsInstance(line)	
		
		if not result.confidentEnough(confidence_level):
			continue

		xh, yh, xc, yc = result.weightData()
		
		if xc < threshold_internal_segment_evidences:
			continue 

		if yc < threshold_split_read_evidences:
			continue

		if result.true_c_vaf not in true_c_vafs: 
			true_c_vafs.append(result.true_c_vaf)
		results.append(result) 

	true_c_vafs.sort()
	return results, true_c_vafs
		

def ln(x):
	"""Natural logarithm. Returns -inf when x is equal to 0 (or below system precision)"""
	if x < sys.float_info.min:
		return float("-inf")
	return math.log(x)

def f(x, mean, std):
		"""Returns probability of observing x given a discretized normal distribution with mean mu and std of sigma, truncated below 0"""
		return (norm.cdf((x + 1.0 - mean) / std) - norm.cdf((x - mean) / std)) / (1.0 - norm.cdf(-1 * mean / std))
	
class ResultsInstance:
		
	def __init__(self, line):
		"""Creates a results instance given the line in the file"""
		self.values 		= line.split('\t') # all values are stored here
		self.chromosome		= self.values[0][3:]
		self.position 		= int(self.values[1])
		self.length 		= int(self.values[2])
		self.true_h_vaf		= self.readFloat(self.values[3])
		self.true_c_vaf 	= self.readFloat(self.values[4])
		self.posterior_h_vaf 	= [float(self.values[5]), float(self.values[6]), float(self.values[7])] 
		self.global_max_exists	= int(self.values[8])
		self.c_vaf0 		= self.readFloat(self.values[9])
		self.maxlogl0 		= self.readFloat(self.values[10])
		self.ci_95_0 		= [self.readFloat(self.values[11]), self.readFloat(self.values[12])] 
		self.ci_90_0 		= [self.readFloat(self.values[13]), self.readFloat(self.values[14])] 
		self.c_vaf1 		= self.readFloat(self.values[9])
		self.maxlogl1 		= self.readFloat(self.values[16])
		self.ci_95_1 		= [self.readFloat(self.values[17]), self.readFloat(self.values[18])] 
		self.ci_90_1 		= [self.readFloat(self.values[19]), self.readFloat(self.values[20])] 
		self.c_vaf2 		= self.readFloat(self.values[9])
		self.maxlogl2 		= self.readFloat(self.values[22])
		self.ci_95_2 		= [self.readFloat(self.values[23]), self.readFloat(self.values[24])] 
		self.ci_90_2 		= [self.readFloat(self.values[25]), self.readFloat(self.values[26])] 
		self.Xh, self.Xh_p	= self.decodeData(self.values[27])
		self.Yh, self.Yh_p	= self.decodeData(self.values[28])
		self.Xc, self.Xc_p	= self.decodeData(self.values[29])
		self.Yc, self.Yc_p	= self.decodeData(self.values[30])
		self.p_value 		= self.readFloat(self.values[31])

	def readFloat(self, x):
		x = x.rstrip()
		x = x.lstrip()
		if x == "None":
			return None
		elif x == "-inf":
			return float("-inf")
		else:
			return float(x)

	def decodeData(self, data_vector):
		"""Decodes the data_vector [ <value_1/align_prob_1> <value_2/align_prob_2> ... <value_n/align_prob_n>]
	   	   Returns a list with the values and a list with the associated alignment probabilities. """
		data_vector 	= data_vector[1:-1] # throw away brackets []
		values 		= [] # list with data values  
		align_probs 	= [] # list with associated alignment probabilities
		for data_point in data_vector.split(' '):
			if data_point == "" or data_point == '[' or data_point == ']':
				continue
			d = data_point.split('/')
			values.append(float(d[0]))
			align_probs.append(float(d[1]))
	
		return values, align_probs
	

	def loglikelihood(self, h_vaf, c_vaf, alpha = 0.0, epsilon0 = 0.0, epsilon1 = 0.0, mu = 112.0, sigma = 15.0):
		"""Returns the loglikelihood. h_vaf and c_vaf denote the VAF of the healthy (h) and cancer (c) cells, respectively."""
		logl = 0.0
		delta = float(self.length)

		# First, we process the alignments of the control sample
		for i in range(len(self.Xh)):
			logl += ln(self.Xh_p[i] * (h_vaf * f(self.Xh[i], mu + delta, sigma) + (1.0 - h_vaf) * f(self.Xh[i], mu, sigma)) + (1.0 - self.Xh_p[i])) 
		for i in range(len(self.Yh)):
			if self.Yh[i] == 1:
				logl += ln(self.Yh_p[i] * (h_vaf*(1.0 - epsilon1) + (1.0 - h_vaf)*epsilon0) + (1.0 - self.Yh_p[i]))
			else:
				logl += ln(self.Yh_p[i] * (h_vaf*epsilon1 + (1.0 - h_vaf)*(1.0 - epsilon0)) + (1.0 - self.Yh_p[i])) 
		# Secondly, we process the alignments of the disease/cancer sample
		for i in range(len(self.Xc)):
			f_delta = f(self.Xc[i], mu + delta, sigma)
			f_0 = f(self.Xc[i], mu, sigma)
			logl += ln(self.Xc_p[i] * (alpha * (h_vaf * f_delta + (1.0 - h_vaf) * f_0) + (1.0 - alpha) * (c_vaf * f_delta + (1.0 - c_vaf) * f_0)) + (1.0 - self.Xc_p[i])) 
		for i in range(len(self.Yc)):
			if self.Yc[i] == 1:
				logl += ln(self.Yc_p[i] * (alpha * (h_vaf*(1.0 - epsilon1) + (1.0 - h_vaf)*epsilon0) + (1.0 - alpha)*(c_vaf*(1.0 - epsilon1) + (1.0 - c_vaf)*epsilon0)) + (1.0 - self.Yc_p[i]))
			else:
				logl += ln(self.Yc_p[i] * (alpha*(h_vaf*epsilon1 + (1.0 - h_vaf)*(1.0 - epsilon0)) + (1.0 - alpha)*((c_vaf*epsilon1 + (1.0 - c_vaf)*(1.0 - epsilon0)))) + (1.0 - self.Yc_p[i])) 
		return logl 
	
	def determineObservedFisherInformation(self, alpha = 0.0, epsilon0 = 0.0, epsilon1 = 0.0, mu = 112.0, sigma = 15.0):
		"""Returns the observed Fisher information. The healthy VAF is set to 0.0."""
		# TODO IMPLEMENT ALPHA!
		self.observed_fisher_information = 0.0
		delta = float(self.length)

		theta_hat = self.c_vaf0 

		for i in range(len(self.Xc)):
			f_present = f(self.Xc[i], mu + delta, sigma)
			f_absent = f(self.Xc[i], mu, sigma)
			pi = self.Xc_p[i]
			self.observed_fisher_information += ((pi * (f_present - f_absent)) / (pi * (theta_hat * f_present + (1.0 - theta_hat) * f_absent) + (1 - pi)))**2
		for i in range(len(self.Yc)):
			f_present = (1.0 - epsilon1)**self.Yc[i] * epsilon1**(1.0 - self.Yc[i])
			f_absent = epsilon0**self.Yc[i] ** (1 - epsilon0)**(1.0 - self.Yc[i])
			pi = self.Yc_p[i]
			self.observed_fisher_information += ((pi * (f_present - f_absent)) / (pi * (theta_hat * f_present + (1.0 - theta_hat) * f_absent) + (1 - pi)))**2

		return self.observed_fisher_information
			
			


	def determinePvalue(self, alpha = 0.0, epsilon0 = 0.0, epsilon1 = 0.0, mu = 112.0, sigma = 15.0):
		self.p_value = self.returnPvalue(alpha = 0.0, epsilon0 = 0.0, epsilon1 = 0.0, mu = 112.0, sigma = 15.0)

	def returnPvalue(self, alpha = 0.0, epsilon0 = 0.0, epsilon1 = 0.0, mu = 112.0, sigma = 15.0):
		h_vaf_mle = 0.0 
		c_vaf_mle = self.c_vaf0

		T = 2.0*(self.loglikelihood(0.0, c_vaf_mle, alpha = 0.0, epsilon0 = 0.0, epsilon1 = 0.0, mu = 112.0, sigma = 15.0) - self.loglikelihood(0.0,0.0, alpha = 0.0, epsilon0 = 0.0, epsilon1 = 0.0, mu = 112.0, sigma = 15.0)) 
		p = 1.0 - chi2.cdf(T, 1) 
		return p 

	def fallsInLengthRange (self, lower = None, upper = None):
		"""Checks whether an indel is in a length range. In case of 'None', that side of the interval is considered unbounded, e.g., [250, None] are all indels of length larger than 250"""
		if lower != None:
			if self.length < lower:
				return False
		if upper != None:
			if self.length > upper:
				return False
		return True
	
	def fallsInCancerVAFRange (self, lower = None, upper = None):
		"""Checks whether the true cancer VAF lies within a certain interval. In case of 'None', that side of the interval is considered unbounded."""
		if lower != None:
			if self.true_c_vaf < lower:
				return False
		if upper != None:
			if self.true_c_vaf > upper:
				return False
		return True

	def confidentEnough(self, confidence_level = 0.0):
		"""Determines whether there is a genotype for the control which is above the confidence level"""
		if max(self.posterior_h_vaf) >= confidence_level:
			return True
		return False

	def numberOfDataPoints(self):
		"""Returns the total number of data points of the different types."""
		return len(self.Xh),  len(self.Yh), len(self.Xc), len(self.Yc)      

	def weightData (self):
		"""Returns the sum of the alignment probabilities for the different types and sources of data."""
		xh, yh, xc, yc = 0.0, 0.0, 0.0, 0.0
		for p in self.Xh_p: xh += p
		for p in self.Yh_p: yh += p
		for p in self.Xc_p: xc += p
		for p in self.Yc_p: yc += p
		return xh, yh, xc, yc 

	def returnNumberSplitReads (self):
		"""Returns the sum of the alignment probabilities of the overlapping alignments that support the presence of the indel."""
		a, b = 0.0, 0.0
		for i in range(len(self.Yh)):
			if self.Yh[i] == 1:
				a += self.Yh_p[i] 
		for i in range(len(self.Yc)):
			if self.Yc[i] == 1:
				b += self.Yc_p[i] 	
		return a, b	
	

	def genotypeControl(self, confidence_level = 0.0):
		"""Returns the genotype (0.0, 0.5, 1.0) of the healthy VAF on the basis of the control sample. Returns 'None' when there is not enough confidence."""
		if not self.confidentEnough(confidence_level = confidence_level):
			return None # not confident enough to make any statement
		if self.posterior_h_vaf[0] >= self.posterior_h_vaf[1] and self.posterior_h_vaf[0] >= self.posterior_h_vaf[2]:
			return 0.0 # homozygous for not having the indel
		if self.posterior_h_vaf[1] >= self.posterior_h_vaf[0] and self.posterior_h_vaf[1] >= self.posterior_h_vaf[2]:
			return 0.5 # heterozygous 
		return 1.0 # homozygous for having the indel
		

	def returnMLEs(self):
		"""Returns the MLEs of the healthy and the cancer cells."""
		h_vaf_mle = 0.0
		c_vaf_mle = self.c_vaf0
		if self.maxlogl1 > self.maxlogl0 and self.maxlogl1 > self.maxlogl2:
			h_vaf_mle = 0.5
			c_vaf_mle = self.c_vaf1
		if self.maxlogl2 > self.maxlogl0 and self.maxlogl2 > self.maxlogl1:
			h_vaf_mle = 1.0
			c_vaf_mle = self.c_vaf2
		return h_vaf_mle, c_vaf_mle

	def returnCancerVAF_MLE(self):
		"""Returns the cancer VAF MLE and the 95% confidence interval."""
		c_vaf_mle = self.c_vaf0
		ci = self.ci_95_0
		if self.maxlogl1 > self.maxlogl0 and self.maxlogl1 > self.maxlogl2:
			c_vaf_mle = self.c_vaf1
			ci = self.ci_95_1
		if self.maxlogl2 > self.maxlogl0 and self.maxlogl2 > self.maxlogl1:
			c_vaf_mle = self.c_vaf2
			ci = self.ci_95_2
		return c_vaf_mle, ci

	def call(self, p_value_threshold = 0.05, confidence_level = 0.0, confidence_interval = 0.95, width_conf_interval = 1.0, use_ci = False):
		"""Returns 'somatic', 'germline' or 'not present'. Returns None if there is not enough certainty about the call"""
		if not self.confidentEnough(confidence_level=confidence_level) or self.global_max_exists == 0:
			return None
	
		h_vaf_mle = 0.0
		c_vaf_mle = self.c_vaf0
		ci 	  = self.ci_95_0
		if confidence_interval is not .95:
			ci = self.ci_90_0 

		if self.maxlogl1 > self.maxlogl0 and self.maxlogl1 > self.maxlogl2:
			h_vaf_mle = 0.5
			c_vaf_mle = self.c_vaf1
			ci 	  = self.ci_95_1
			if confidence_interval is not .95:
				ci = self.ci_90_1 
		if self.maxlogl2 > self.maxlogl0 and self.maxlogl2 > self.maxlogl1:
			h_vaf_mle = 1.0
			c_vaf_mle = self.c_vaf2
			ci 	  = self.ci_95_2
			if confidence_interval is not .95:
				ci = self.ci_90_2 

		
		if ci[1] - ci[0] > width_conf_interval:
			return None

		if h_vaf_mle > 0.0:
			return 'germline'

		if not use_ci:
			if self.p_value < p_value_threshold: 
				return 'somatic'
			else:
				return 'not present'
		else:
			if confidence_interval == 0.95:
				if self.ci_95_0[0] == 0.0 and self.ci_95_0[1] == 1.0:
					return None # not confident enough  

				if self.ci_95_0[0] == 0.0:
					return 'not present'
				else:
					return 'somatc'
			if confidence_interval == 0.90:
				if self.ci_90_0[0] == 0.0 and self.ci_90_0[1] == 1.0:
					return None # not confident enough 

				if self.ci_90_0[0] == 0.0:
					return 'not present'
				else:
					return 'somatic'
	

	def obtainTruth(self):
		"""Returns 'somatic', 'germline' or 'not present'"""
		if self.true_h_vaf > 0.0:
			return 'germline'
		if self.true_c_vaf > 0.0:
			return 'somatic'
		return 'not present'

	def isSomatic(self):
		"""Returns True is the genetic variant is somatic."""
		if self.true_h_vaf == 0.0 and self.true_c_vaf > 0:
			return True
		return False

	def isGermline(self):
		"""Returns True is the genetic variant is germline."""
		if self.true_h_vaf > 0.0 and self.true_c_vaf == 0:
			return True
		return False

	def isNotPresent(self):
		"""Returns True is the genetic variant is not present, i.e., not germline and not somatic"""
		if self.true_h_vaf == 0.0 and self.true_c_vaf == 0:
			return True
		return False

	def print(self):
		print(self.values)

