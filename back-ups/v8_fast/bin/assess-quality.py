#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
from sets import Set
import os
import sys
import math

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__))[:-3] + 'python')
from Indel import * 
from ContingencyTable import * 

__author__ = "Louis Dijkstra"

usage = """%prog <annotated-vcf-file>

	<annotated-vcf-file> 	A VCF file processes the somatic mutation caller.

Assesses the quality of the classifications made by the somatic mutation caller.

The analysis is done for various length ranges, which can be changed inside the file;
look for LENGTH_RANGES. 
"""

LENGTH_RANGES 		= [[None, None], [None,9], [10,29], [30,49], [50,99], [100, None], [250, None]]
CANCER_VAF_RANGES 	= [[None, None], [0.0, 0.25], [.25, .50], [.50, .75], [.75, None]]
GAMMA_LEVELS 		= [10.0**(-x) for x in range(10)] 

def returnAutosome (autosome):
	"""Returns an integer denoting the chromosome."""
	if len(autosome) > 3 and autosome[:3] == 'chr':
		autosome = autosome[3:]
	if autosome == 'X' or autosome == 'x':
		return 23
	if autosome == 'Y' or autosome == 'y':
		return 24	
	return int(autosome)

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
	
def returnMinimumDifference (list1, list2):
	"""Returns minimal difference between all possible pairs of elements contained in list1 and list2."""
	min_diff = float("inf")
	for e1 in list1:
		for e2 in list2:
			if min_diff > abs(e1 - e2):
				min_diff = abs(e1 - e2)
	return min_diff





class TestMSE:
	"""Assesses the quality of the c_vaf estimate."""
	def __init__(self):
		self.se_del = [0.0 for i in range(len(LENGTH_RANGES))] 
		self.n_del = [0.0 for i in range(len(LENGTH_RANGES))] 
		self.se_ins = [0.0 for i in range(len(LENGTH_RANGES))] 
		self.n_ins = [0.0 for i in range(len(LENGTH_RANGES))] 

	def update(self, t_deletion, t_length, t_index, t_c_vaf, c_est_c_vaf):
		for l in range(len(LENGTH_RANGES)):
			length_range = LENGTH_RANGES[l]
			for i in range(len(t_deletion)):
				if not fallsInInterval(t_length[i], length_range):
					continue 
				if len(t_index[i]) == 0:
					continue
				est = None
				for j in range(len(t_index[i])):
					est = c_est_c_vaf[t_index[i][j]]
					if est != None:
						break ; 
				if est == None:
					continue
				if t_deletion[i] == 1:
					self.se_del[l] += (t_c_vaf[i] - est)**2.0 
					self.n_del[l] += 1 
				else:
					self.se_ins[l] += (t_c_vaf[i] - est)**2.0 
					self.n_ins[l] += 1 

	def printResults(self):
		print('\n\nQUALTIY OF CANCER VAF ESTIMATE (MEAN SQUARED ERROR)\n')
		print('DELETION')
		print('length range\tmean squared error')
		print('------------\t------------------')
		for l in range(len(LENGTH_RANGES)):
			print(intervalToStr(LENGTH_RANGES[l]), '\t', end = '')
			if self.n_del[l] == 0:
				print('-')
			else: 
				print(self.se_del[l] / self.n_del[l])
		print('------------\t------------------\n')
		print('INSERTION')
		print('length range\tmean squared error')
		print('------------\t------------------')
		for l in range(len(LENGTH_RANGES)):
			print(intervalToStr(LENGTH_RANGES[l]), '\t', end = '')
			if self.n_ins[l] == 0:
				print('-')
			else: 
				print(self.se_ins[l] / self.n_ins[l])
		print('------------\t------------------\n')

class TestGenotyping: 
	"""Assesses the quality of the h_vaf estimate."""
	def __init__(self):
		self.tables_del = [Table3x3() for i in range(len(LENGTH_RANGES))]  
		self.tables_ins = [Table3x3() for i in range(len(LENGTH_RANGES))]  
		
	def update(self, t_deletion, t_length, t_index, t_h_vaf, c_est_h_vaf):
		for l in range(len(LENGTH_RANGES)):
			length_range = LENGTH_RANGES[l]
			for i in range(len(t_deletion)):
				if not fallsInInterval(t_length[i], length_range):
					continue 
				if len(t_index[i]) == 0:
					continue
				row = 0
				column = 0 
				est = c_est_h_vaf[t_index[i][0]]
				if (t_h_vaf[i] == 0.5): column = 1
				if (t_h_vaf[i] == 0.0): column = 2
				if (est == 0.5): row = 1
				if (est == 0.0): row = 2
				if t_deletion[i] == 1:
					self.tables_del[l].add(row, column)
				else:
					self.tables_ins[l].add(row, column)
		
	def printResults(self):
		print('\n\nQUALTIY OF HEALTHY VAF ESTIMATE (GENOTYPING)\n')
		print('DELETION')
		print('length range\taccuracy') 
		print('------------\t--------')
		for l in range(len(LENGTH_RANGES)):
			print(intervalToStr(LENGTH_RANGES[l]), '\t', end = '')
			print(self.tables_del[l].returnPercentageCorrect())
		print('------------\t--------\n')
		print('INSERTION')
		print('length range\taccuracy') 
		print('------------\t--------')
		for l in range(len(LENGTH_RANGES)):
			print(intervalToStr(LENGTH_RANGES[l]), '\t', end = '')
			print(self.tables_ins[l].returnPercentageCorrect())
		print('------------\t--------\n')


class TestBinary: 
	
	def __init__(self, label, title):
		self.classification = label
		self.title = title

		self.d_n = [[0 for c in range(len(CANCER_VAF_RANGES))] for l in range(len(LENGTH_RANGES))] 
		self.d_tp = [[[0 for g in range(len(GAMMA_LEVELS))] for c in range(len(CANCER_VAF_RANGES))] for l in range(len(LENGTH_RANGES))] 
		self.d_call = [[[0 for g in range(len(GAMMA_LEVELS))] for c in range(len(CANCER_VAF_RANGES))] for l in range(len(LENGTH_RANGES))] 
		self.d_call_tp = [[[0 for g in range(len(GAMMA_LEVELS))] for c in range(len(CANCER_VAF_RANGES))] for l in range(len(LENGTH_RANGES))] 
		self.i_n = [[0 for c in range(len(CANCER_VAF_RANGES))] for l in range(len(LENGTH_RANGES))] 
		self.i_tp = [[[0 for g in range(len(GAMMA_LEVELS))] for c in range(len(CANCER_VAF_RANGES))] for l in range(len(LENGTH_RANGES))] 
		self.i_call = [[[0 for g in range(len(GAMMA_LEVELS))] for c in range(len(CANCER_VAF_RANGES))] for l in range(len(LENGTH_RANGES))] 
		self.i_call_tp = [[[0 for g in range(len(GAMMA_LEVELS))] for c in range(len(CANCER_VAF_RANGES))] for l in range(len(LENGTH_RANGES))] 
	
	
	def update(self, t_deletion, t_length, t_index, t_c_vaf, t_class, c_deletion, c_length, c_index, c_est_c_vaf, c_class, c_posterior_prob):
		for l in range(len(LENGTH_RANGES)):
			length_range = LENGTH_RANGES[l] 
			for c in range(len(CANCER_VAF_RANGES)):
				cancer_vaf_range = CANCER_VAF_RANGES[c]
				# walk through the true indels
				for i in range(len(t_deletion)):
					if not fallsInInterval(t_length[i], length_range):
						continue
					if not fallsInInterval(t_c_vaf[i], cancer_vaf_range):
						continue
					if t_class[i] != self.classification:
						continue
					if t_deletion[i] == 1:	self.d_n[l][c] += 1
					else:			self.i_n[l][c] += 1
					classification = None
					post_prob = None
					for j in range(len(t_index[i])):
						classification = c_class[t_index[i][j]]
						post_prob = c_posterior_prob[t_index[i][j]]
						if classification != None:
							break 
					if classification != self.classification:
						continue
					for g in range(len(GAMMA_LEVELS)):
						threshold = 1.0 - GAMMA_LEVELS[g]
						if post_prob >= threshold:
							if t_deletion[i] == 1:	self.d_tp[l][c][g] += 1
							else:			self.i_tp[l][c][g] += 1
				# walk through the calls
				for i in range(len(c_deletion)):
					if not fallsInInterval(c_length[i], length_range):
						continue
					if not fallsInInterval(c_est_c_vaf[i], cancer_vaf_range):
						continue
					if c_class[i] != self.classification:
						continue
					for g in range(len(GAMMA_LEVELS)):
						threshold = 1.0 - GAMMA_LEVELS[g] 
						if c_posterior_prob[i] >= threshold: 
							if c_deletion[i] == 1:	self.d_call[l][c][g] += 1
							else:			self.i_call[l][c][g] += 1
							classification = None
							for j in range(len(c_index[i])):
								classification = t_class[c_index[i][j]]
								if classification != None:
									break 
							if classification == self.classification:
								if c_deletion[i] == 1:	self.d_call_tp[l][c][g] += 1
								else:			self.i_call_tp[l][c][g] += 1

	def printResults(self):			
		print("\n\n", self.title, "\n")
		print("DELETIONS")
		print("length range\tcancer vaf range\tgamma\t# true som.\t# som. calls\trecall\tprecision")
		print("------------\t----------------\t-----\t-----------\t------------\t------\t---------")
		for l in range(len(LENGTH_RANGES)):
			for c in range(len(CANCER_VAF_RANGES)):
				for g in range(len(GAMMA_LEVELS)):
					print(intervalToStr(LENGTH_RANGES[l]), '\t', end = '')
					print(intervalToStr(CANCER_VAF_RANGES[c]), '\t\t', end = '')
					print(GAMMA_LEVELS[g], '\t\t', end = '')
					print(self.d_n[l][c], '\t', self.d_call[l][c][g], '\t\t', end = '')
					recall_string, precision_string = '-.----', '-.----'	
					if self.d_n[l][c] != 0:
						recall_string = "{:.4f}".format(float(self.d_tp[l][c][g]) / float(self.d_n[l][c]))
					if self.d_call[l][c][g] != 0:
						precision_string = "{:.4f}".format(float(self.d_call_tp[l][c][g]) / float(self.d_call[l][c][g]))
					print(recall_string, '\t', precision_string)
				print("------------\t----------------\t-----\t-----------\t------------\t------\t---------")
		
		print("\nDELETIONS -- NO BOUNDS ON THE CANCER VAF")
		print("length range\tgamma\t# true som.\t# som. calls\trecall\tprecision")	
		print("------------\t-----\t-----------\t------------\t------\t---------")
		for l in range(len(LENGTH_RANGES)):
			for g in range(len(GAMMA_LEVELS)):
				print(intervalToStr(LENGTH_RANGES[l]), '\t', end = '')
				print(GAMMA_LEVELS[g], '\t\t', end = '')
				print(self.d_n[l][0], '\t', self.d_call[l][0][g], '\t\t', end = '')
				recall_string, precision_string = '-.----', '-.----'	
				if self.d_n[l][0] != 0:
					recall_string = "{:.4f}".format(float(self.d_tp[l][0][g]) / float(self.d_n[l][0]))
				if self.d_call[l][0][g] != 0:
					precision_string = "{:.4f}".format(float(self.d_call_tp[l][0][g]) / float(self.d_call[l][0][g]))
				print(recall_string, '\t', precision_string)	
			print("------------\t-----\t-----------\t------------\t------\t---------")
		
		print("\n\nINSERTIONS")
		print("length range\tcancer vaf range\tgamma\t# true som.\t# som. calls\trecall\tprecision")
		print("------------\t----------------\t-----\t-----------\t------------\t------\t---------")
		for l in range(len(LENGTH_RANGES)):
			for c in range(len(CANCER_VAF_RANGES)):
				for g in range(len(GAMMA_LEVELS)):
					print(intervalToStr(LENGTH_RANGES[l]), '\t', end = '')
					print(intervalToStr(CANCER_VAF_RANGES[c]), '\t\t', end = '')
					print(GAMMA_LEVELS[g], '\t', end = '')
					print(self.i_n[l][c], '\t\t', self.i_call[l][c][g], '\t\t', end = '')
					recall_string, precision_string = '-.----', '-.----'	
					if self.i_n[l][c] != 0:
						recall_string = "{:.4f}".format(float(self.i_tp[l][c][g]) / float(self.i_n[l][c]))
					if self.i_call[l][c][g] != 0:
						precision_string = "{:.4f}".format(float(self.i_call_tp[l][c][g]) / float(self.i_call[l][c][g]))
					print(recall_string, '\t', precision_string)
				print("------------\t----------------\t-----\t-----------\t------------\t------\t---------")
	
		print("\nINSERTIONS -- NO BOUNDS ON THE CANCER VAF")
		print("length range\tgamma\t# true som.\t# som. calls\trecall\tprecision")	
		print("------------\t-----\t-----------\t------------\t------\t---------")
		for l in range(len(LENGTH_RANGES)):
			for g in range(len(GAMMA_LEVELS)):
				print(intervalToStr(LENGTH_RANGES[l]), '\t', end = '')
				print(GAMMA_LEVELS[g], '\t\t', end = '')
				print(self.i_n[l][0], '\t', self.i_call[l][0][g], '\t\t', end = '')
				recall_string, precision_string = '-.----', '-.----'	
				if self.i_n[l][0] != 0:
					recall_string = "{:.4f}".format(float(self.i_tp[l][0][g]) / float(self.i_n[l][0]))
				if self.i_call[l][0][g] != 0:
					precision_string = "{:.4f}".format(float(self.i_call_tp[l][0][g]) / float(self.i_call[l][0][g]))
				print(recall_string, '\t', precision_string)	
			print("------------\t-----\t-----------\t------------\t------\t---------")
		print("")


class TestSomatic(TestBinary):
	
	def __init__(self):
		TestBinary.__init__(self, 'SOMATIC', 'SOMATIC/NOT SOMATIC') 
		
class TestGermline(TestBinary):
	
	def __init__(self):
		TestBinary.__init__(self, 'GERMLINE', 'GERMLINE/NOT GERMLINE') 

class TestAbsent(TestBinary):
	
	def __init__(self):
		TestBinary.__init__(self, 'ABSENT', 'ABSENT/PRESENT') 

class TestCaller: 
	""""Asesses the quality of the caller"""
	def __init__(self):
		self.del_n_true_positives 	= [0 for i in range(len(LENGTH_RANGES))] ; # indels called correctly
		self.del_n 			= [0 for i in range(len(LENGTH_RANGES))] ; # total number of indels
		self.del_n_false_positives 	= [0 for i in range(len(LENGTH_RANGES))] ; # indels called falsely
		self.del_n_calls		= [0 for i in range(len(LENGTH_RANGES))] ; # total number of calls made
		self.ins_n_true_positives 	= [0 for i in range(len(LENGTH_RANGES))] ; # indels called correctly
		self.ins_n 			= [0 for i in range(len(LENGTH_RANGES))] ; # total number of indels
		self.ins_n_false_positives 	= [0 for i in range(len(LENGTH_RANGES))] ; # indels called falsely
		self.ins_n_calls		= [0 for i in range(len(LENGTH_RANGES))] ; # total number of calls made

	def update(self, t_deletion, t_length, t_index, c_deletion, c_length, c_index):
		for l in range(len(LENGTH_RANGES)):
			length_range = LENGTH_RANGES[l]			
			for i in range(len(t_deletion)):
				if not fallsInInterval(t_length[i], length_range):
					continue 
				if t_deletion[i] == 1:
					self.del_n[l] += 1
					if len(t_index[i]) != 0:
						self.del_n_true_positives[l] += 1
				else: 
					self.ins_n[l] += 1
					if len(t_index[i]) != 0:
						self.ins_n_true_positives[l] += 1
		for l in range(len(LENGTH_RANGES)):
			length_range = LENGTH_RANGES[l]			
			for i in range(len(c_deletion)):
				if not fallsInInterval(c_length[i], length_range):
					continue 
				if c_deletion[i] == 1:
					self.del_n_calls[l] += 1
					if len(c_index[i]) == 0:
						self.del_n_false_positives[l] += 1
				else: 
					self.ins_n_calls[l] += 1
					if len(c_index[i]) == 0:
						self.ins_n_false_positives[l] += 1
			 
	def summarizeTruth(self):
		"""Summarizes the truth found in the true VCF file."""
		print('\n\nSUMMARY\n')
		print('length range\t# true deletions\t# true insertions\t# true indels')
		print('------------\t----------------\t-----------------\t-------------')
		for l in range(len(LENGTH_RANGES)):
			print(intervalToStr(LENGTH_RANGES[l]), '\t', self.del_n[l], '\t\t\t', self.ins_n[l], '\t\t\t', self.del_n[l] + self.ins_n[l]) 
		print('------------\t----------------\t-----------------\t-------------')

	def printResults(self):		
		"""Prints the results to the command line."""
		print('\n\nPERFORMANCE OF THE CALLER\n')
		
		print('DELETIONS')
		print('length range\t# deletions\trecall\tprecision') 
		print('------------\t-----------\t------\t---------')
		for l in range(len(LENGTH_RANGES)):
			print(intervalToStr(LENGTH_RANGES[l]), '\t', end = '')
			print(self.del_n[l], end = '\t\t')
			recall_string, precision_string = '-.----', '-.----'		
			if self.del_n[l] != 0:
				recall_string = "{:.4f}".format(float(self.del_n_true_positives[l]) / float(self.del_n[l]))
			if self.del_n_calls[l] != 0:
				precision_string = "{:.4f}".format(float(self.del_n_calls[l] - self.del_n_false_positives[l]) / float(self.del_n_calls[l]))		
			print(recall_string, '\t', precision_string) 
		print('------------\t-----------\t------\t---------\n\n')	
			

		print('INSERTIONS')
		print('length range\t# insertions\trecall\tprecision') 
		print('------------\t------------\t------\t---------')
		for l in range(len(LENGTH_RANGES)):
			print(intervalToStr(LENGTH_RANGES[l]), '\t', end = '')
			print(self.ins_n[l], end = '\t\t')
			recall_string, precision_string = '-.----', '-.----'		
			if self.ins_n[l] != 0:
				recall_string = "{:.4f}".format(float(self.ins_n_true_positives[l]) / float(self.ins_n[l]))
			if self.ins_n_calls[l] != 0:
				precision_string = "{:.4f}".format(float(self.ins_n_calls[l] - self.ins_n_false_positives[l]) / float(self.ins_n_calls[l]))		
			print(recall_string, '\t', precision_string) 
		print('------------\t------------\t------\t---------\n')	

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

	true_class = 'GERMLINE'
	if true_h_vaf == 0.0:
		true_class = 'SOMATIC'

	return true_h_vaf, true_c_vaf, true_class 			


def returnTitleTable (length_range, confidence_level, cancer_vaf_range, confidence_interval):
	"""Returns an appropriate description for a table given the parameters used."""
	if confidence_interval == None: 
		return 'length: ' + returnIntervalString(length_range) + " | confidence level: " + str(confidence_level) + " | cancer vaf range: " + returnIntervalString(cancer_vaf_range) 
	else: 
		return 'length: ' + returnIntervalString(length_range) + " | confidence level: " + str(confidence_level) + " | cancer vaf range: " + returnIntervalString(cancer_vaf_range) + " | confidence interval used: " + str(confidence_interval * 100) + "%"

	
def readFloat(x):
	x = x.strip()
	if x == "NONE":	return None
	if x == "None":	return None
	if x == "none": return None
	if x == "unknown": return None
	if x == "UNKNOWN": return None
	if x == "-inf":	return float("-inf")
	return float(x)

def returnClassification(vcf_record):
	map_h_vaf 	= readFloat(vcf_record.INFO['MAP_HEALTHY_VAF'][0])
	map_c_vaf 	= readFloat(vcf_record.INFO['MAP_CANCER_VAF'][0])
	call 		= vcf_record.INFO['CALL'][0].strip()
	post_prob 	= readFloat(vcf_record.INFO['POSTERIOR_PROB'][0])
	return map_h_vaf, map_c_vaf, call, post_prob 


def main():

	parser = OptionParser(usage=usage)
	parser.add_option("--even", action="store_true", dest="even_chromosomes", default=False,
						help="Summarize the results for all even chromosomes, while normally only the uneven ones are used. ONLY USE IN THE END! This construction is used to avoid overfitting.")
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
	parser.add_option("-t", action="store", dest="truth_summary_file", default=None, 
				  		help="The location of file that contains the ground truth. (Default = data/truth/truth.summary")
	parser.add_option("-x", action="store", dest="autosome", default=None, type=int, 
				  		help="Processes just one autosome.")
	(options, args) = parser.parse_args()
	
	if (len(args)!=1):
		parser.print_help()
		return 1

	# get truth file
	if options.truth_summary_file == None:
		options.truth_summary_file = os.path.abspath(os.path.dirname(__file__))[:-3] + 'data/truth/truth.summary'

	if options.high_precision:
		options.length_threshold 	= 10
		options.distance_threshold 	= 20

	list_of_autosomes = range(1,23,2) # only the uneven autosomes
	if options.even_chromosomes: 
		list_of_autosomes = range(2,23,2)
	if options.all_chromosomes:
		list_of_autosomes = range(1,23)
	if options.autosome != None:
		list_of_autosomes = [options.autosome]
	
	print('INPUT:\tfile with ground truth @', os.path.abspath(options.truth_summary_file), '\n\tcalls in VCF file @', os.path.abspath(args[0]), '\n')

	truth_file = open(options.truth_summary_file, 'r') # true indels 
	call_vcf_reader = vcf.Reader(open(os.path.abspath(args[0]))) # called indels

	test_caller = TestCaller()
	test_genotyping = TestGenotyping()
	test_mse = TestMSE()
	test_absent = TestAbsent()
	test_germline = TestGermline()
	test_somatic = TestSomatic()

	# work through the true VCF file:
	for autosome in list_of_autosomes: 
			# allocate memory 
		t_deletion 		= [] 
		t_centerpoints 		= [] 
		t_length		= []
		t_h_vaf			= []
		t_c_vaf 		= [] 
		t_class 		= [] 
		t_index		 	= [] 	

		c_deletion		= []
		c_centerpoints		= []
		c_length		= []
		c_est_h_vaf		= []
		c_est_c_vaf		= [] 
		c_class			= [] 
		c_posterior_prob 	= [] 
		c_index 		= [] 
		
		print("Processing autosome %d"%autosome) 
		
		try:
			values = truth_file.next().split('\t')
		
			# read on until you reach the variant on the wanted autosome
			while (int(values[0]) != autosome):
				try:
					values = truth_file.next().split('\t')
				except:
					break 
		
			print("Reading in VCF records from %s"%options.truth_summary_file)
			while (int(values[0]) == autosome):
				t_deletion.append(int(values[1]))
				t_length.append(int(values[2]))
				t_h_vaf.append(readFloat(values[3]))
				t_c_vaf.append(readFloat(values[4]))
				t_class.append(values[5].strip())
				centerpoints = [int(values[6])]
				if len(values) == 8:
					centerpoints.append(int(values[7]))
				t_centerpoints.append(centerpoints)
				t_index.append([]) 
				try:
					values = truth_file.next().split('\t')
				except:
					break 
		except:
			pass

		

		c_vcf_record = call_vcf_reader.next() 
		
		# read on until you reach the variant on the wanted autosome
		while (returnAutosome(c_vcf_record.CHROM) != autosome):			
			c_vcf_record = call_vcf_reader.next() 
			
		print("Reading in VCF records from %s"%args[0])
		while (returnAutosome(c_vcf_record.CHROM) == autosome):
			is_deletion = isDeletion(c_vcf_record)
			is_insertion = isInsertion(c_vcf_record)
			if (not is_deletion and not is_insertion): # not an indel
				c_vcf_record = call_vcf_reader.next() 
			else:
				length = int(returnIndelLength(c_vcf_record))
				c_length.append(length)
				centerpoints = [] 
				if is_deletion:	
					c_deletion.append(1) 
					if length % 2 == 0: 
						centerpoints = [c_vcf_record.POS + (length/2), c_vcf_record.POS + length / 2 + 1]
					else: 
						centerpoints = [c_vcf_record.POS + int(length / 2)] 
				else: 		
					c_deletion.append(0)
					centerpoints = [c_vcf_record.POS, c_vcf_record.POS + 1]
				c_centerpoints.append(centerpoints) 	
				call_h_vaf, call_c_vaf, call_class, p = returnClassification(c_vcf_record) 	
				c_est_h_vaf.append(call_h_vaf)
				c_est_c_vaf.append(call_c_vaf)
				c_class.append(call_class) 
				c_posterior_prob.append(p)
				c_index.append([]) 
				c_vcf_record = call_vcf_reader.next() 

		print("Sorting the VCF records from %s"%args[0])
		c_centerpoints, c_deletion, c_length, c_est_h_vaf, c_est_c_vaf, c_class, c_posterior_prob, c_index = (list(t) for t in zip(*sorted(zip(c_centerpoints, c_deletion, c_length, c_est_h_vaf, c_est_c_vaf, c_class, c_posterior_prob,  c_index))))

		print("Match the records from both files for autosome %d"%autosome) 
		low, up = 0,0 
		for i in range(len(t_centerpoints)):
			# update low and up
			while (low < len(c_centerpoints) and min(t_centerpoints[i]) - options.distance_threshold > min(c_centerpoints[low])):
				low += 1

			up = low 
			while (up < len(c_centerpoints) and max(t_centerpoints[i]) + options.distance_threshold >= max(c_centerpoints[up])):
				up += 1 

			for j in range(low, up):
				if t_deletion[i] != c_deletion[j]:
					continue  
				if returnMinimumDifference(t_centerpoints[i], c_centerpoints[j]) <= options.distance_threshold and abs(t_length[i] - c_length[j]) <= options.length_threshold:
					t_index[i].append(j)
					c_index[j].append(i)
				
					
		print("Done matching the records from both files for autosome %d"%autosome) 

		print("Update results for autosome %d"%autosome)
		test_caller.update(t_deletion, t_length, t_index, c_deletion, c_length, c_index)
		test_genotyping.update(t_deletion, t_length, t_index, t_h_vaf, c_est_h_vaf)	
		test_mse.update(t_deletion, t_length, t_index, t_c_vaf, c_est_c_vaf)
		test_somatic.update(t_deletion, t_length, t_index, t_c_vaf, t_class, c_deletion, c_length, c_index, c_est_c_vaf, c_class, c_posterior_prob)
		test_germline.update(t_deletion, t_length, t_index, t_c_vaf, t_class, c_deletion, c_length, c_index, c_est_c_vaf, c_class, c_posterior_prob)
		test_absent.update(t_deletion, t_length, t_index, t_c_vaf, t_class, c_deletion, c_length, c_index, c_est_c_vaf, c_class, c_posterior_prob)


	test_caller.summarizeTruth()
	test_caller.printResults() 
	test_genotyping.printResults() 
	test_mse.printResults()
	test_absent.printResults()
	test_germline.printResults()
	test_somatic.printResults()
	
if __name__ == '__main__':
	sys.exit(main())
