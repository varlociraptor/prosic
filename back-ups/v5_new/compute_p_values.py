#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
import os
import sys
import math
from scipy.stats import norm 
from scipy.stats import chi2 

__author__ = "Louis Dijkstra"

usage = """%prog <results-filename>

Computes the p-values and stores those in the last column.

	<results-filename> 	File with the somatic mutation results. 
	
NOTE: 	At the moment we only process cases where the cancer VAF is assumed
		continuous. In addition, only autosomes are taken into account.
"""	

def ln(x):
	"""Natural logarithm. Returns -inf when x is equal to 0 (or below system precision)"""
	if x < sys.float_info.min:
		return float("-inf")
	return math.log(x)

def f(x, mean, std):
		"""Returns probability of observing x given a discretized normal distribution with mean mu and std of sigma, truncated below 0"""
		return (norm.cdf((x + 1.0 - mean) / std) - norm.cdf((x - mean) / std)) / (1.0 - norm.cdf(-1 * mean / std))

class OldResultsInstance:
		
	def __init__(self, line):
		"""Creates a results instance given the line in the file"""
		self.values 		= line.split('\t') # all values are stored here
		self.chromosome		= self.values[0][3:]
		self.position 		= int(self.values[1])
		self.length 		= int(self.values[2])
		self.true_h_vaf		= float(self.values[3])
		self.true_c_vaf 	= float(self.values[4])
		self.posterior_h_vaf 	= [float(self.values[5]), float(self.values[6]), float(self.values[7])] 
		self.global_max_exists	= int(self.values[8])
		self.c_vaf0 		= self.readFloat(self.values[9])
		self.maxlogl0 		= float(self.values[10])
		self.ci_95_0 		= [float(self.values[11]), float(self.values[12])] 
		self.ci_90_0 		= [float(self.values[13]), float(self.values[14])] 
		self.c_vaf1 		= self.readFloat(self.values[9])
		self.maxlogl1 		= float(self.values[16])
		self.ci_95_1 		= [float(self.values[17]), float(self.values[18])] 
		self.ci_90_1 		= [float(self.values[19]), float(self.values[20])] 
		self.c_vaf2 		= self.readFloat(self.values[9])
		self.maxlogl2 		= float(self.values[22])
		self.ci_95_2 		= [float(self.values[23]), float(self.values[24])] 
		self.ci_90_2 		= [float(self.values[25]), float(self.values[26])] 
		self.Xh, self.Xh_p	= self.decodeData(self.values[27])
		self.Yh, self.Yh_p	= self.decodeData(self.values[28])
		self.Xc, self.Xc_p	= self.decodeData(self.values[29])
		self.Yc, self.Yc_p	= self.decodeData(self.values[30])

	def readFloat(self, x):
		if x == "None" or x == "None ":
			return None
		elif x == "-inf" or x == "-inf ":
			return float("-inf")
		elif x == None:
			return None
		else:
			return float(x)

	def decodeData(self, data_vector):
		"""Decodes the data_vector [ <value_1/align_prob_1> <value_2/align_prob_2> ... <value_n/align_prob_n>]
	   	   Returns a list with the values and a list with the associated alignment probabilities. """
		data_vector = data_vector[1:-1] # throw away brackets []
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
	
	def determinePvalue(self, alpha = 0.0, epsilon0 = 0.0, epsilon1 = 0.0, mu = 112.0, sigma = 15.0):
		self.p_value = self.returnPvalue(alpha = 0.0, epsilon0 = 0.0, epsilon1 = 0.0, mu = 112.0, sigma = 15.0)

	def returnPvalue(self, alpha = 0.0, epsilon0 = 0.0, epsilon1 = 0.0, mu = 112.0, sigma = 15.0):
		h_vaf_mle = 0.0 
		c_vaf_mle = self.c_vaf0

		T = 2.0*(self.loglikelihood(0.0, c_vaf_mle, alpha = 0.0, epsilon0 = 0.0, epsilon1 = 0.0, mu = 112.0, sigma = 15.0) - self.loglikelihood(0.0,0.0, alpha = 0.0, epsilon0 = 0.0, epsilon1 = 0.0, mu = 112.0, sigma = 15.0)) 
		p = 1.0 - chi2.cdf(T, 1) 
		return p 



def main():

	parser = OptionParser(usage=usage)
	(options, args) = parser.parse_args()
	
	if (len(args)!=1):
		parser.print_help()
		return 1

	results_filename = os.path.abspath(args[0])
	results_file = open(results_filename, 'r')

	results = [] # list with all the results

	# walk through all the results in the file
	for line in results_file: 
		result = OldResultsInstance(line)
		result.determinePvalue()
		for i in range(len(result.values) - 1):
			print(result.values[i], '\t', end = '')
		print(result.values[-1].rstrip(), '\t', end = '')
		print(result.p_value)
		
				
if __name__ == '__main__':
	sys.exit(main())
