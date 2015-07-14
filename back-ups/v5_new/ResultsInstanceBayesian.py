#!/usr/bin/env python
from __future__ import print_function, division
from scipy.stats import norm 
from scipy.stats import chi2 

__author__ = "Louis Dijkstra"

"""
	ResultsInstanceBayesian.py contains the functionailty to process intermediate results when a pure Bayesian approach is used. 
"""

def obtainResults(results_file, all_chromosomes = False, confidence_level = 0.0, threshold_internal_segment_evidences = 0.0, threshold_split_read_evidences = 0.0):
	"""Returns all the results from a results file from chromosome 1, 2 and 3. When the boolean all_chromosomes = True, then all results are returned."""
	results 	= [] # list with all the results
	true_c_vafs 	= [] # a list of all the true different cancer VAFs present in the file

	truth_vcf_file_used = False
	if 'truth' in results_file.name:
		truth_vcf_file_used = True
	
	# walk through all the results in the file
	for line in results_file: 

		if len(line.split('\t')) < 26:
			continue

		chrom = line[3:5]
		if truth_vcf_file_used:
			chrom = line[:2]

		if chrom == 'X' or chrom == 'X ' or chrom == 'Y' or chrom == 'Y ' or chrom == '': # sex chromosomes are not taken into account
				continue
 
		if not all_chromosomes:			
			if int(chrom) > 3:
				continue 

		result = ResultsInstanceBayesian(line)	
		
		if not result.confidentEnoughAboutGenotype(confidence_level):
			continue

		if not result.confidentEnoughAboutCalling(confidence_level):
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
	
class ResultsInstanceBayesian:
		
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
		self.p_somatic_mutation = self.readFloat(self.values[9])
		self.p_germline 	= self.readFloat(self.values[10])
		self.p_not_present	= self.readFloat(self.values[11])
		self.ci			= [self.readFloat(self.values[12]), self.readFloat(self.values[13])]
		self.mle_h_vaf		= self.readFloat(self.values[14])
		self.mle_c_vaf		= self.readFloat(self.values[15])
		self.c_vaf0 		= self.readFloat(self.values[16])
		self.maxlogl0 		= self.readFloat(self.values[17])
		self.c_vaf1 		= self.readFloat(self.values[18])
		self.maxlogl1 		= self.readFloat(self.values[19])
		self.c_vaf2 		= self.readFloat(self.values[20])
		self.maxlogl2 		= self.readFloat(self.values[21])
		self.Xh, self.Xh_p	= self.decodeData(self.values[22])
		self.Yh, self.Yh_p	= self.decodeData(self.values[23])
		self.Xc, self.Xc_p	= self.decodeData(self.values[24])
		self.Yc, self.Yc_p	= self.decodeData(self.values[25])

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

	def confidentEnoughAboutGenotype(self, confidence_level = 0.0):
		"""Determines whether there is a genotype for the control which is above the confidence level"""
		if max(self.posterior_h_vaf) >= confidence_level:
			return True
		return False

	def confidentEnoughAboutCalling(self, confidence_level = 0.0):
		"""Determines whether there is a hypothesis (somatic/germline/not present) which is above the confidence level"""
		if self.p_somatic_mutation >= confidence_level or self.p_germline >= confidence_level or self.p_not_present >= confidence_level:
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
		if not self.confidentEnoughAboutGenotype(confidence_level = confidence_level):
			return None # not confident enough to make any statement
		if self.posterior_h_vaf[0] >= self.posterior_h_vaf[1] and self.posterior_h_vaf[0] >= self.posterior_h_vaf[2]:
			return 0.0 # homozygous for not having the indel
		if self.posterior_h_vaf[1] >= self.posterior_h_vaf[0] and self.posterior_h_vaf[1] >= self.posterior_h_vaf[2]:
			return 0.5 # heterozygous 
		return 1.0 # homozygous for having the indel


	def call(self, confidence_level = 0.0):
		"""Returns 'somatic', 'germline' or 'not present'. Returns None if there is not enough certainty about the call"""
		if not self.confidentEnoughAboutCalling(confidence_level=confidence_level) or self.global_max_exists == 0:
			return None
	
		if self.p_somatic_mutation > self.p_germline and self.p_somatic_mutation > self.p_not_present:
			return 'somatic'
		if self.p_germline > self.p_somatic_mutation and self.p_germline > self.p_not_present:
			return 'germline'
		return 'not present'

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

