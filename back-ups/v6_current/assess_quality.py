#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
from sets import Set
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
		self.true_positives 	= [] # indels called correctly
		self.false_negatives 	= [] # indels that were not called
		#self.n_called_deletions = 0 # total number of deletions called

		self.correct_calls = []
		self.correct_calls_pos = [] 

		self.list_of_autosomes = list_of_autosomes

	#def setNCalledDeletions(self, n_called_deletions):
	#	self.n_called_deletions = n_called_deletions

	def addCorrectCall(self, vcf_record):
		if not vcf_record.POS in self.correct_calls_pos:
			self.correct_calls_pos.append(vcf_record.POS)
			self.correct_calls.append(vcf_record)

	def addTruePositive(self, call):
		true_vcf_record 	= call[0]
		similar_vcf_record 	= call[1]
		self.true_positives.append(true_vcf_record)
		self.addCorrectCall(similar_vcf_record)

	def addFalseNegative(self, true_vcf_record):
		self.false_negatives.append(true_vcf_record)

	def returnNumberIndels(self, list_vcf_records, length_range, cancer_vaf_range):
		"""Returns the number of deletions and insertions in <list_vcf_records> that fall
		   within in the length and cancer VAF range."""
		n_true_deletions, n_true_insertions = 0, 0 
		for true_vcf_record in list_vcf_records: 
			if not returnAutosome(true_vcf_record.CHROM) in self.list_of_autosomes:
				continue
			indel = Indel(true_vcf_record)
			if not fallsInInterval(indel.length, length_range):
				continue
			true_h_vaf, true_c_vaf, true_class = returnTruth(true_vcf_record)
			if not fallsInInterval(true_c_vaf, cancer_vaf_range):
				continue
			if isDeletion(true_vcf_record): n_true_deletions 	+= 1
			if isInsertion(true_vcf_record): n_true_insertions 	+= 1
		return n_true_deletions, n_true_insertions

	def returnNumberTrueIndels(self, length_range, cancer_vaf_range):
		"""Returns the number of true deletions and insertions."""
		n_true_deletions, n_true_insertions = self.returnNumberIndels(self.true_positives, length_range, cancer_vaf_range)
		n_d, n_i = self.returnNumberIndels(self.false_negatives, length_range, cancer_vaf_range)
		return n_true_deletions + n_d, n_true_insertions + n_i
		


	def returnRecallPrecision(self, length_range):
		"""Determines the recall and precision for deletions from the true VCF file that lie 
		   within a given length range."""
		pass
		
	def summarizeTruth(self):
		"""Summarizes the truth found in the true VCF file."""
		print('\n\nSUMMARY\n')
		print('length range\tcancer VAF range\t# true deletions\t# true insertions\t# true indels')
		print('------------\t----------------\t----------------\t-----------------\t-------------')
		for l in range(len(LENGTH_RANGES)):
			for c in range(len(CANCER_VAF_RANGES)):
				n_deletions, n_insertions = self.returnNumberTrueIndels(LENGTH_RANGES[l], CANCER_VAF_RANGES[c])
				print(intervalToStr(LENGTH_RANGES[l]), '\t', intervalToStr(CANCER_VAF_RANGES[c]), '\t\t', n_deletions, '\t\t', n_insertions, '\t\t', n_deletions + n_insertions)
			print('------------\t----------------\t----------------\t-----------------\t-------------')

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

	def store(self, true_vcf_record, similar_vcf_record, true_h_vaf, map_h_vaf):
		row, column = 0, 0
		if true_h_vaf == 0.5: 	column = 1
		if true_h_vaf == 1.0: 	column = 2
		if map_h_vaf == 0.5: 	row = 1
		if map_h_vaf == 1.0: 	row = 2

		is_deletion = isDeletion(true_vcf_record)
		length = returnIndelLength(true_vcf_record)
		
		self.addResult([is_deletion, length, true_h_vaf], row, column)
	
	def printResults(self, print_tables = False):
		"""Prints the results to standard output"""

		print('\nMAIN TEST - Calling Somatic/Germline/Not Present\n')
		print('length range\tconfidence level\t% correct\t# deletions\ttotal in that range')
		print('------------\t----------------\t---------\t-----------\t-------------------')
		for length_range in LENGTH_RANGES:
			for confidence_level in CONFIDENCE_LEVELS:  
				percentage_correctly_called, n_test_cases_used, total, table = testCallingPerformance(results, length_range, confidence_level, CANCER_VAF_RANGE)	
				perc = "--.----"
				if percentage_correctly_called is not None:				
					perc = "{:.4f}".format(percentage_correctly_called * 100)
				print(returnIntervalString(length_range), '\t', confidence_level, '\t\t\t', perc, '\t', n_test_cases_used, '\t\t', total)
				print('------------\t----------------\t---------\t-----------\t-------------------')
	
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
		self.deletion 		= []
		self.length		= []
		self.true_c_vaf		= []  
		self.squared_errors  	= [] 

	def store(self, true_vcf_record, true_c_vaf, map_c_vaf):
		"""Stores a result."""
		if isDeletion(true_vcf_record):
			self.deletion.append(True)
		else:
			self.deletion.append(False)
		self.length.append(returnIndelLength(true_vcf_record))
		self.true_c_vaf.append(true_c_vaf)
		self.squared_errors.append((true_c_vaf - map_c_vaf)**2.0)

	def returnMSE(self, length_range, cancer_vaf_range):
		"""Computes the MSE for deletions, insertions and all indels for which the 
		   true length and the true cancer VAF fall within the given intervals."""

		# initialize 
		n_del, n_ins, sse_del, sse_ins = 0.0, 0.0, 0.0, 0.0  # sse - sum squared error

		# walk through the stored data 
		for i in range(len(self.deletion)): 
			if not fallsInInterval(length[i], length_range):
				continue
			if not fallsInInterval(true_c_vaf[i], cancer_vaf_range):
				continue
			if deletion[i]:
				sse_del += self.squared_errors[i]
				n_del 	+= 1.0
			else:
				sse_ins += self.squared_errors[i]
				n_ins 	+= 1.0
		return n_del, n_ins, sse_del/n_del, sse_ins/n_ins, (sse_del + sse_ins)/(n_del + n_ins), 
		

	def printResults(self):
		"""Prints the results to standard output"""
		print('\n**************************************************************************************************************')
		print('Mean squared error of the MAP estimate of the cancer VAF')
		print('**************************************************************************************************************\n')
		print('length range\tcancer VAF range\t# deletions\tMSE (deletions)\t# insertions\tMSE (insertions)\t# indel\tMSE (total)')
		print('------------\t----------------\t-----------\t---------------\t------------\t----------------\t-------\t-----------')
		for l in range(len(LENGTH_RANGES)):
			for c in range(len(CANCER_VAF_RANGES)):
				n_d, n_i, mse_d, mse_i, mse = self.returnMSE(LENGTH_RANGES[l], CANCER_VAF_RANGES[c])
				print(intervalToStr(LENGTH_RANGES[l]), '\t', intervalToStr(CANCER_VAF_RANGES[c]), '\t\t', n_d, '\t\t', mse_d, '\t\t', '\t\t', n_i, '\t\t', mse_i, '\t\t', n_i + n_d, '\t\t', mse)
			print('------------\t----------------\t-----------\t---------------\t------------\t----------------\t-------\t-----------')


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
						help="Summarize the results for all chromosomes, while normally only the uneven ones are used. ONLY USE IN THE END! This construction is used to avoid overfitting.")
	parser.add_option("--table", action="store_true", dest="print_tables", default=False,
						help="Prints the 3x3 tables for genotyping and somatic/germline/not present classification.")
	parser.add_option("--high", action="store_true", dest="high_precision", default=False,
						help="Uses high precision when matching detected deletions with the deletions in the true VCF file, i.e., the centerpoints must be less or exactly 20 bp away and the length must not differ more than 10 bp. This option overwrites the options '-d' and '-l' (!).")
	parser.add_option("-d", action="store", dest="distance_threshold", default=100, type=int,
				  		help="Distance between centerpoints of deletions threshold. The distance must be smaller than this in order for the deletions to be deemed similar (in addition, the lengths must be similar, see option '-l'). (Default = 100 bp)")
	parser.add_option("-l", action="store", dest="length_threshold", default=100, type=int,
				  		help="Lengths of deletions threshold. The length difference must be smaller than this in order for the deletions to be deemed similar (in addition, the placement of the centerpoints must be similar, see option '-d'). (Default = 100 bp)")
	parser.add_option("-t", action="store", dest="true_vcf_file", default="data/truth.indels.sorted.vcf.gz", 
				  		help="The location of the VCF file with that contains the true deletions. (Default = data/truth.indels.sorted.vcf.gz")
	(options, args) = parser.parse_args()
	
	if (len(args)!=1):
		parser.print_help()
		return 1

	if options.high_precision:
		options.length_threshold 	= 10
		options.distance_threshold 	= 20

	# TODO change back!
	list_of_autosomes = [1]
	#list_of_autosomes = range(1,23,2) # only the uneven autosomes
	if options.all_chromosomes:
		list_of_autosomes = range(1,23)
	
	print('INPUT:\ttrue VCF file @', os.path.abspath(options.true_vcf_file), '\n\tcalls in VCF file @', os.path.abspath(args[0]), '\n')

	vcf_reader = vcf.Reader(open(options.true_vcf_file)) # true indels
	# indels that are called: 
	vcf_matcher = VCFMatcher(os.path.abspath(args[0]), length_threshold = options.length_threshold, centerpoint_dist_threshold = options.distance_threshold)

	# the tests that are run: 
	caller_test 	= CallerTest(list_of_autosomes) 
	genotype_test 	= GenotypeTest()
	main_test 	= MainTest()
	mse_results 	= MSETest() 

	# determines the total number of deletions called 
	# caller_test.setNCalledDeletions(vcf_matcher.retrieveNumberOfDeletions(list_of_autosomes = list_of_autosomes)) 

	# walk through all true deletions n = 0
	n = 0 
	for true_vcf_record in vcf_reader: 
		n += 1
		if n % 1000 == 0:
			print('Processed %d VCF records'%n, file=sys.stderr)

		if not returnAutosome(true_vcf_record.CHROM) in list_of_autosomes:
			continue
		if not isIndel(true_vcf_record):
			continue		
		
		# retrieve a similar vcf record from the calls made
		similar_vcf_record = vcf_matcher.retrieveMatchingRecord(true_vcf_record, search_range=10000)
		
		if similar_vcf_record == None: # no similar indel...
			caller_test.addFalseNegative(true_vcf_record) 
			continue

		caller_test.addTruePositive([true_vcf_record, similar_vcf_record]) 

		# obtain truth 
		true_h_vaf, true_c_vaf, true_class = returnTruth(true_vcf_record) 		

		# obtain classification
		map_h_vaf, map_c_vaf, called_class = returnClassification(similar_vcf_record) 

		# add result to the main test
		#main_test.store([true_vcf_record, similar_vcf_record], true_class, called_class)

		# add result to the genotyping test 
		#genotype_test.store([true_vcf_record, similar_vcf_record], true_h_vaf, map_h_vaf)
		
		# compute MSE
		mse_results.store(true_vcf_record, true_c_vaf, map_c_vaf) 

	caller_test.summarizeTruth() # prints info on the number of indels of certain length ranges etc.

	#caller_test.printResults()
	#genotype_test.printResults(print_tables = options.print_tables)
	#main_test.printResults(print_tables = options.print_tables)
	mse_results.printResults()
		
		
	
if __name__ == '__main__':
	sys.exit(main())
