#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
from scipy.stats import norm 
from scipy.stats import chi2 
import os
import sys
import math

from Variant.Indel import * 
from VCFMatcher import * 
from ContingencyTable import *

__author__ = "Louis Dijkstra"

usage = """%prog <annotated-vcf-file>

Assesses the quality of the classifications made by the somatic mutation caller in 
'call_somatic_mutations.py'

The analysis is done for various length ranges, which can be changed inside the file;
look for LENGTH_RANGES. 

	<annotated-vcf-file> 	A VCF file processes by 'call_somatic_mutations.py'
"""

LENGTH_RANGES 		= [[None, None], [10,29], [30,49], [50,69], [70,99], [100, 249], [250, None]]
CANCER_VAF_RANGES 	= [[None, None], [0.0, 0.20], [.20, .40], [.40, .60], [.60, .80], [.80, 1.0], [.20, None], [.40, None], [.50, None], [.60, None]]
CONFIDENCE_LEVELS 	= [0.0, .50, .90, .95, .99]

def fallsInInterval(value, interval):
	"""Returns True when value lies within in the interval."""
	if interval[0] != None:
		if value < interval[0]:
			return False
	if interval[1] != None:
		if value > interval[1]:
			return False
	return True

def intervalToStr(interval):
	if interval[0] is None and interval[1] is not None:
		return "at most " + str(interval[1])
	if interval[0] is not None and interval[1] is None:
		return "at least " + str(interval[0])
	if interval[0] is None and interval[1] is None:
		return "unbounded"
	return "[" + str(interval[0]) + ','  + str(interval[1]) + "]"

def returnFancyIntervalString(interval):
	"""Returns the interval in the form of a string. Suitable for plotting."""
	if interval[0] is None:
		if interval[1] is None: # unbounded
			return r".-."
		else:	# <= interval[1]
			return r"$\leq$" + str(interval[1])
	elif interval[1] is None: # >= interval[0]
			return r"$\geq$" + str(interval[0])
	return str(interval[0]) + '-'  + str(interval[1]) 
	
class Test3x3: 
	"""Class for storing the results of a test with three classes."""
	def __init__(self, labels = ['0.0', '0.5', '1.0']):
		self.labels = labels
		self.results = [[[] for i in range(3)] for j in range(3)]
	
	def addResult(self, obj, row, column):
		self.results[row][column].append(obj)

class CallerTest:
	"""Assesses the quality of the caller; whether it called the right deletions."""
	def __init__(self, list_of_autosomes):
		# allocate memory
		self.true_positives = [] # deletions called correctly
		self.false_negatives = [] # deletions that were not called
		self.n_called_deletions = 0 # total number of deletions called

		self.list_of_autosomes = list_of_autosomes

	def setNCalledDeletions(self, n_called_deletions):
		self.n_called_deletions = n_called_deletions

	def addTruePositive(self, true_deletion):
		self.true_positives.append(true_deletion)

	def addFalseNegative(self, true_deletion):
		self.false_negatives.append(true_deletion)

	def returnNumberTrueDeletions(self, length_range, cancer_vaf_range):
		n_true_deletions = 0
		for deletion in self.true_positives:
			if not returnAutosome(deletion.vcf_record.CHROM) in self.list_of_autosomes:
				continue 
			if not fallsInInterval(deletion.length, length_range):
				continue
			true_h_vaf, true_c_vaf, true_class = returnTruth(deletion.vcf_record)
			if not fallsInInterval(true_c_vaf, cancer_vaf_range):
				continue
			n_true_deletions += 1 
		for deletion in self.false_negatives:
			if not returnAutosome(deletion.vcf_record.CHROM) in self.list_of_autosomes:
				continue 
			if not fallsInInterval(deletion.length, length_range):
				continue
			true_h_vaf, true_c_vaf, true_class = returnTruth(deletion.vcf_record)
			if not fallsInInterval(true_c_vaf, cancer_vaf_range):
				continue
			n_true_deletions += 1
		return n_true_deletions 

	def returnRecallPrecision(self, length_range):
		"""Determines the recall and precision for deletions from the true VCF file that lie 
		   within a given length range."""
		
		
	def summarizeTruth(self):
		"""Summarizes the truth found in the true VCF file."""
		print('SUMMARY\n')
		print('length range\tcancer VAF range\t# true deletions')
		print('------------\t----------------\t----------------')
		for l in range(len(LENGTH_RANGES)):
			for c in range(len(CANCER_VAF_RANGES)):
				print(intervalToStr(LENGTH_RANGES[l]), '\t', intervalToStr(CANCER_VAF_RANGES[c]), '\t\t', self.returnNumberTrueDeletions(LENGTH_RANGES[l], CANCER_VAF_RANGES[c]))
			print('------------\t----------------\t----------------')

	def printResults(self):
		"""Prints the results to standard output"""
		print('RECALL/PRECISION FOR THE CALLER IN QUESTION\n')
		print('length range\trecall\tprecision')
		print('------------\t------\t---------')
		for length_range in LENGTH_RANGES:
			

recall_string, precision_string = '-.----', '-.----'				
			if recall is not None:
				recall_string = "{:.4f}".format(recall)
			if precision is not None:
				precision_string = "{:.4f}".format(precision)
		

class GenotypeTest(Test3x3):
	def __init__(self):
		Test3x3.__init__(self, labels = ['0.0', '0.5', '1.0'])	

	def store(self, true_deletion, called_deletion, true_h_vaf, map_h_vaf):
		row, column = 0, 0
		if true_h_vaf == 0.5: 	column = 1
		if true_h_vaf == 1.0: 	column = 2
		if map_h_vaf == 0.5: 	row = 1
		if map_h_vaf == 1.0: 	row = 2
		self.addResult([true_deletion, called_deletion], row, column)
	
	def printResults(self, print_tables = False):
		"""Prints the results to standard output"""
		pass

class MainTest(Test3x3):
	def __init__(self):
		Test3x3.__init__(self, labels = ['somatic', 'germline', 'not present'])
			
	def store(self, true_deletion, called_deletion, true_h_vaf, map_h_vaf):
		row, column = 0, 0
		if true_class == 'germline': 		column = 1
		if true_class == 'not present': 	column = 2
		if called_class == 'germline': 		row = 1
		if called_class == 'not present': 	row = 2
		self.addResult([true_deletion, called_deletion], row, column)

	def printResults(self, print_tables = False):
		"""Prints the results to standard output"""
		pass

class MSETest:

	def __init__(self):
		self.objects 	= []
		self.values 	= [] 

	def store(self, true_deletion, called_deletion, true_c_vaf, map_c_vaf):
		self.objects.append([true_deletion, called_deletion])
		self.values.append((true_c_vaf - map_c_vaf)**2.0)

	def printResults(self):
		"""Prints the results to standard output"""
		pass

def determineVAF(gt_nums):
	if gt_nums == '1|1':
		return 1.0 
	elif gt_nums == None:
		return 0.0 
	return 0.5
			
def returnTruth(true_vcf_record):
	"""Returns the true VAFs of the healthy and cancer cells and the classification (somatic/germline)"""
	# There is one control VAF and VAFs for the four cancer populations
	true_h_vaf, som1_vaf, som2_vaf, som3_vaf, som4_vaf = 0.0, 0.0, 0.0, 0.0, 0.0 
	for call in true_vcf_record.samples:
		if call.sample=='Control': 	true_h_vaf 	= determineVAF(call.gt_nums)
		if call.sample=='Som1':		som1_vaf 	= determineVAF(call.gt_nums)	
		if call.sample=='Som2':		som2_vaf 	= determineVAF(call.gt_nums)	
		if call.sample=='Som3':		som3_vaf 	= determineVAF(call.gt_nums)
		if call.sample=='Som4':		som4_vaf 	= determineVAF(call.gt_nums)
	true_c_vaf = (1/3.0)*som1_vaf + (1/3.0)*som2_vaf + (1/4.0)*som3_vaf + (1/12.0)*som4_vaf 

	true_class = 'germline'
	if true_h_vaf == 0.0:
		true_class = 'somatic'

	return true_h_vaf, true_c_vaf, true_class 			


def returnTitleTable (length_range, confidence_level, cancer_vaf_range, confidence_interval):
	"""Returns an appropriate description for a table given the parameters used."""
	if confidence_interval == None: 
		return 'length: ' + returnIntervalString(length_range) + " | confidence level: " + str(confidence_level) + " | cancer vaf range: " + returnIntervalString(cancer_vaf_range) 
	else: 
		return 'length: ' + returnIntervalString(length_range) + " | confidence level: " + str(confidence_level) + " | cancer vaf range: " + returnIntervalString(cancer_vaf_range) + " | confidence interval used: " + str(confidence_interval * 100) + "%"

	
def readFloat(x):
	x = x.rstrip()
	x = x.lstrip()
	if x == "None":	return None
	if x == "-inf":	return float("-inf")
	return float(x)

def returnClassification(similar_vcf_record):
	map_h_vaf 	= readFloat(similar_vcf_record.INFO['MAP_HEALTHY_VAF'])
	map_c_vaf 	= readFloat(similar_vcf_record.INFO['MAP_CANCER_VAF'])
	call 		= similar_vcf_record.INFO['CALL']


def main():

	parser = OptionParser(usage=usage)
	parser.add_option("--final", action="store_true", dest="all_chromosomes", default=False,
						help="Summarize the results for all chromosomes, while normally only the first three are used. ONLY USE IN THE END! This construction is used to avoid overfitting.")
	parser.add_option("--table", action="store_true", dest="print_tables", default=False,
						help="Prints the 3x3 tables for genotyping and somatic/germline/not present classification.")
	parser.add_option("--high", action="store_true", dest="high_precision", default=False,
						help="Uses high precision when matching detected deletions with the deletions in the true VCF file, i.e., the centerpoints must be less or exactly 20 bp away and the length must not differ more than 10 bp. This option overwrites the options '-d' and '-l' (!).")
	parser.add_option("-d", action="store", dest="distance_threshold", default=100, type=int,
				  		help="Distance between centerpoints of deletions threshold. The distance must be smaller than this in order for the deletions to be deemed similar (in addition, the lengths must be similar, see option '-l'). (Default = 100 bp)")
	parser.add_option("-l", action="store", dest="length_threshold", default=100, type=int,
				  		help="Lengths of deletions threshold. The length difference must be smaller than this in order for the deletions to be deemed similar (in addition, the placement of the centerpoints must be similar, see option '-d'). (Default = 100 bp)")
	parser.add_option("-t", action="store", dest="true_vcf_file", default="truth.vcf.gz", 
				  		help="The location of the VCF file with that contains the true deletions. (Default = ~/Projects/SomaticMutationCalling/Data/truth.vcf.gz)")
	(options, args) = parser.parse_args()
	
	if (len(args)!=1):
		parser.print_help()
		return 1

	if options.high_precision:
		options.length_threshold = 10
		options.distance_threshold = 20

	list_of_autosomes = range(1,4) # only the first three chromosomes
	if options.all_chromosomes:
		list_of_autosomes = range(1,23)
	
	print('INPUT:\ttrue VCF file @', os.path.abspath(options.true_vcf_file), '\n\tcalls in VCF file @', os.path.abspath(args[0]), '\n')

	vcf_reader = vcf.Reader(open(options.true_vcf_file)) # true deletions 
	# deletions that are called: 
	vcf_matcher = VCFMatcher(os.path.abspath(args[0]), length_threshold = options.length_threshold, centerpoint_dist_threshold = options.distance_threshold)

	caller_test 	= CallerTest(list_of_autosomes) 
	genotype_test 	= GenotypeTest()
	main_test 	= MainTest()
	mse_results 	= MSETest() 

	# determines the total number of deletions called 
	# caller_test.setNCalledDeletions(vcf_matcher.retrieveNumberOfDeletions(list_of_autosomes = list_of_autosomes)) 

	# walk through all true deletions 
	for true_vcf_record in vcf_reader: 
		if not returnAutosome(true_vcf_record.CHROM) in list_of_autosomes:
			continue

		true_deletion = Deletion(true_vcf_record)

		# obtain the list of similar deletions in the input VCF file
		similar_vcf_records = vcf_matcher.retrieveMatchingDeletions(true_deletion, search_range = 10000)	

		if len(similar_vcf_records) == 0: # no deletions similar...  
			caller_test.addFalseNegative(true_deletion) 
			continue

		caller_test.addTruePositive(true_deletion) 
		#similar_vcf_record = similar_vcf_records[0]
		#if len(similar_vcf_records) > 1: # in case there are more similar deletions
		#	similar_vcf_record = selectMostSimilarDeletion(true_deletion, similar_vcf_records) # TODO
		
		# obtain truth 
		#true_h_vaf, true_c_vaf, true_class = returnTruth(true_vcf_record) 
		
		# obtain classification
		#similar_deletion = Deletion(similar_vcf_record)
		#map_h_vaf, map_c_vaf, called_class = returnClassification(similar_vcf_record) 

		# add result to the main test
		#main_test.store(true_deletion, similar_deletion, true_class, called_class)

		# add result to the genotyping test 
		#genotype_test.store(true_deletion, similar_deletion, true_h_vaf, map_h_vaf)

		# compute MSE
		#mse_results.store(true_deletion, similar_deletion, true_c_vaf, map_c_vaf) 
	
	caller_test.summarizeTruth() # prints info on the number of deletions of certain length ranges etc.

	#caller_test.printResults()
	#genotype_test.printResults(print_tables = options.print_tables)
	#main_test.printResults(print_tables = options.print_tables)
	#mse_results.printResults()
		
		
	
if __name__ == '__main__':
	sys.exit(main())
