#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
from scipy.stats import norm 
from scipy.stats import chi2 
import os
import sys
import math

from ResultsInstance import *
from Indel import * 

__author__ = "Louis Dijkstra"

usage = """%prog <results-filename>

	<results-filename> 	File with the somatic mutation results. 
	
NOTE: 	At the moment we only process cases where the cancer VAF is assumed
		continuous. In addition, only autosomes are taken into account.
"""

LENGTH_RANGES 		= [[None, None], [10,29], [30,49], [50,69], [70,99], [100, 249], [250, None]]
CANCER_VAF_RANGES 	= [[None, None], [0.0, 0.20], [.20, .40], [.40, .60], [.60, .80], [.80, 1.0], [.20, None], [.40, None], [.50, None], [.60, None]]
CONFIDENCE_LEVELS 	= [0.0, .50, .90, .95, .99]
CANCER_VAF_RANGE 	= [None, None]

VCF_TRUTH_FILENAME = "/data1/structvar/vaf-experiments/vcf/cancerclones.vcf" # VCF file that contains the ground truth
VCF_TRUTH_READER = vcf.Reader(open(VCF_TRUTH_FILENAME)) # used for reading the VCF file 

list_true_deletions = [[] for i in range(24)] # 24 empty lists for the 22 autosomes and the two sex-chromosomes X and Y   


def somatic(coverage, true_vcf_record):
	"""Returns True when the true VCF record represents a somatic mutation, False otherwise."""
	true_h_vaf, true_c_vaf = returnTrueVAFs (coverage, true_vcf_record)
	if true_h_vaf == 0 and true_c_vaf > 0:
		return True
	return False

def determineVAF (gt_nums):
	if gt_nums == '1|1':
		return 1.0 
	elif gt_nums == None:
		return 0.0 
	return 0.5

def returnTrueVAFs (coverage, true_vcf_record):
	"""Returns the true VAFs of the healthy and cancer cells."""
	# There is one control VAF and VAFs for the four cancer populations
	true_h_vaf, som1_vaf, som2_vaf, som3_vaf, som4_vaf = 0.0, 0.0, 0.0, 0.0, 0.0 
	for call in true_vcf_record.samples:
		if call.sample=='Control': 	true_h_vaf 	= determineVAF(call.gt_nums)
		if call.sample=='Som1':		som1_vaf 	= determineVAF(call.gt_nums)	
		if call.sample=='Som2':		som2_vaf 	= determineVAF(call.gt_nums)	
		if call.sample=='Som3':		som3_vaf 	= determineVAF(call.gt_nums)
		if call.sample=='Som4':		som4_vaf 	= determineVAF(call.gt_nums)
		
	if coverage == 40:
		return true_h_vaf, (1/3.0)*som1_vaf + (1/3.0)*som2_vaf + (1/6.0)*som3_vaf + (1/6.0)*som4_vaf 
	else:
		return true_h_vaf, (1/3.0)*som1_vaf + (1/3.0)*som2_vaf + (1/4.0)*som3_vaf + (1/12.0)*som4_vaf 

def returnIndex (chromosome):
	if len(chromosome) > 3 and chromosome[:3] == 'chr':
		chromosome = chromosome[3:]
	if chromosome == 'X' or chromosome == 'x':
		return 22
	if chromosome == 'Y' or chromosome == 'y':
		return 23	
	return int(chromosome)

def obtainListOfTrueDeletions (min_length): 
	"""Fills the global variable list_true_deletions with the true deletions from the true VCF file"""
	for true_vcf_record in VCF_TRUTH_READER:
		if not isDeletion(true_vcf_record): # not a deletion
			continue 
		if returnIndelLength(true_vcf_record) < min_length: # too small
			continue 
		list_true_deletions[returnIndex(true_vcf_record.CHROM)].append(true_vcf_record) 

def deletionsSimilar(deletion, result, difference_length_threshold = 100, difference_centerpoint_threshold = 100):
	"""Returns true when deletion and result are similar, otherwise false."""
	# determine centerpoints
	centerpoints_other_del = [] # centerpoints of the other deletion
	if result.length % 2 == 0: 
		centerpoints_other_del = [result.position + 1 + result.length / 2 - 1, result.position + 1 + result.length / 2]
	else: 
		centerpoints_other_del = [result.position + 1 + result.length / 2]
	if returnMinimumDifference(deletion.centerpoints, centerpoints_other_del) > difference_centerpoint_threshold:
		return False
	if abs(deletion.length - result.length) > difference_length_threshold:
		return False
	return True

def returnRecallPrecision (table2x2): 
	"""Returns the recall and precision given a 2x2 contingency table. Returns 'None' when not defined"""
	recall, precision = None, None
	column_total = table2x2[0][0] + table2x2[1][0]
	row_total = table2x2[0][0] + table2x2[0][1]
	if column_total != 0:
		recall = table2x2[0][0] / float(column_total)
	if row_total != 0:
		precision = table2x2[0][0] / float(row_total)
	return recall, precision

def printTable2x2(table, title = "", labels = ["Somatic", "Not somatic"]):
	"""Prints a 2x2 table to standard output"""
	row_total, column_total, total = [0,0], [0,0], 0
	# compute the marignals and the total
	for r in [0,1]:
		for c in [0,1]:		
			row_total[r] 	+= table[r][c]
			column_total[c] += table[r][c]
			total 		+= table[r][c]
	recall, precision = returnRecallPrecision(table)
	print('------------------------' + title + '------------------------')
	print('\t\t\t\t\ttruth')
	print('\t\t|\t', labels[0], '\t\t|\t', labels[1], '\t|\ttotal')
	print('--------------------------------------------------------------------------------------------')
	print(labels[0], '\t\t|\t', table[0][0], '\t\t|\t', table[0][1], '\t\t|\t', row_total[0])
	print(labels[1], '\t\t|\t', table[1][0], '\t\t|\t', table[1][1], '\t\t|\t', row_total[1])	
	print('--------------------------------------------------------------------------------------------')
	print('total\t\t|\t', column_total[0], '\t\t|\t', column_total[1], '\t\t|\t', total)
	print('\nrecall:\t\t', recall)
	print('precision:\t', precision)	
	print('-------------------------------------------------------------------\n\n')
	

def printTable3x3(table, title = "", labels = ["0.0", "0.5", "1.0"]):
	"""Prints a 3x3 table to standard output"""
	row_total, column_total, total = [0,0,0], [0,0,0], 0
	# compute the marignals and the total
	for r in [0,1,2]:
		for c in [0,1,2]:		
			row_total[r] 	+= table[r][c]
			column_total[c] += table[r][c]
			total 		+= table[r][c]

	n_correct = 0 # number of correctly classified item
	for i in [0,1,2]: 
		n_correct += table[i][i] # on the diagonal

	print('------------------------' + title + '------------------------')
	print('\t\t\t\t\t\ttruth')
	print('\t\t\t|\t', labels[0], '\t|\t', labels[1], '\t|\t', labels[2], '\t|\ttotal')
	print('\t--------------------------------------------------------------------------------------------')
	print('\t\t', labels[0], '\t|\t',table[0][0], '\t|\t', table[0][1],'\t|\t', table[0][2],'\t|\t', row_total[0])
	print('assign.\t\t', labels[1], '\t|\t',table[1][0], '\t|\t', table[1][1],'\t|\t', table[1][2],'\t|\t', row_total[1])
	print('\t\t', labels[2], '\t|\t',table[2][0], '\t|\t', table[2][1],'\t|\t', table[2][2],'\t|\t', row_total[2])
	print('\t--------------------------------------------------------------------------------------------')
	print('\t\ttotal\t|\t', column_total[0], '\t|\t', column_total[1], '\t|\t', column_total[2], '\t|\t', total)
	if total != 0:
		print('\ncorrectly classified:\t', n_correct, ' /', n_correct / float(total) * 100.0, '%\n')
	else: 
		print('\n')


def returnIntervalString (interval):
	if interval[0] is None and interval[1] is not None:
		return "at most " + str(interval[1])
	elif interval[0] is not None and interval[1] is None:
		return "at least " + str(interval[0])
	elif interval[0] is None and interval[1] is None:
		return "unbounded"
	else:
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

def returnTitleTable (length_range, confidence_level, cancer_vaf_range, confidence_interval):
	"""Returns an appropriate description for a table given the parameters used."""
	if confidence_interval == None: 
		return 'length: ' + returnIntervalString(length_range) + " | confidence level: " + str(confidence_level) + " | cancer vaf range: " + returnIntervalString(cancer_vaf_range) 
	else: 
		return 'length: ' + returnIntervalString(length_range) + " | confidence level: " + str(confidence_level) + " | cancer vaf range: " + returnIntervalString(cancer_vaf_range) + " | confidence interval used: " + str(confidence_interval * 100) + "%"

def testCallingPerformance(results, length_range, confidence_level, cancer_vaf_range):
	"""Returns the performance for calling an indel somatic/germline or not present."""
	n_test_cases_used, total = 0, 0
	table = [[0,0,0],[0,0,0],[0,0,0]]
	for result in results:
		if not result.fallsInLengthRange(lower=length_range[0], upper=length_range[1]):
			continue
		#if not result.fallsInCancerVAFRange(lower=cancer_vaf_range[0], upper=cancer_vaf_range[1]):
		#	continue
		total += 1 # item falls into the predefined class (might not be confident enough!)

		call = result.call(confidence_level = confidence_level)
		if call is None: # enough confidence about the genotype of the healthy cells
			continue
			

		n_test_cases_used += 1
		truth = result.obtainTruth()

		#print(truth, '\t', call, '\t', 0.01 / float(len(results)), '\t', result.p_value)

		row, column = 0, 0 
		if call == 'germline': 		row = 1
		if call == 'not present': 	row = 2
		if truth == 'germline': 	column = 1
		if truth == 'not present': 	column = 2
		table[row][column] += 1
	
	if n_test_cases_used == 0:
		return None, 0, total, table
	
	return (table[0][0] + table[1][1] + table[2][2]) / float(n_test_cases_used), n_test_cases_used, total, table	

def testSomaticMutationCallingPerformance(results, length_range, confidence_level, cancer_vaf_range):
	n_test_cases_used, total = 0, 0
	table = [[0,0],[0,0],[0,0]]
	for result in results:
		if not result.fallsInLengthRange(lower=length_range[0], upper=length_range[1]):
			continue
		if not result.fallsInCancerVAFRange(lower=cancer_vaf_range[0], upper=cancer_vaf_range[1]):
			continue
		total += 1 # item falls into the predefined class (might not be confident enough!)

		call = result.call(confidence_level = confidence_level)
		if call is None: # enough confidence about the genotype of the healthy cells
			continue
	
		if call == 'germline' or call == 'not present':
			call = 'not somatic'
	
		n_test_cases_used += 1
		truth = result.obtainTruth()
		if truth == 'germline' or truth == 'not present': 
			truth = 'not somatic'

		row, column = 0, 0 
		if call == 'not somatic': 	row = 1
		if truth == 'not somatic': 	column = 1
	
		table[row][column] += 1

	if n_test_cases_used == 0:
		return None, None, 0, total
	
	recall, precision = returnRecallPrecision (table)
	return recall, precision, n_test_cases_used, total	


def testPresenceDetectionPerformance(results, length_range, confidence_level, cancer_vaf_range):
	n_test_cases_used, total = 0, 0
	table = [[0,0],[0,0],[0,0]]
	for result in results:
		if not result.fallsInLengthRange(lower=length_range[0], upper=length_range[1]):
			continue
		if not result.fallsInCancerVAFRange(lower=cancer_vaf_range[0], upper=cancer_vaf_range[1]):
			continue
		total += 1 # item falls into the predefined class (might not be confident enough!)

		call = result.call(confidence_level = confidence_level)
		if call is None: # enough confidence about the genotype of the healthy cells
			continue
	
		if call == 'germline' or call == 'somatic':
			call = 'present'
	
		n_test_cases_used += 1
		truth = result.obtainTruth()
		if truth == 'germline' or truth == 'somatic': 
			truth = 'present'

		row, column = 0, 0 
		if call == 'not present': 	row = 1
		if truth == 'not present': 	column = 1
	
		table[row][column] += 1

	if n_test_cases_used == 0:
		return None, None, 0, total
	
	recall, precision = returnRecallPrecision (table)
	return recall, precision, n_test_cases_used, total	
	
def testGenotpyingPerformance(results, length_range, confidence_level):
	n_test_cases_used, total = 0, 0
	table = [[0,0,0],[0,0,0],[0,0,0]]
	for result in results:
		if not result.fallsInLengthRange(lower=length_range[0], upper=length_range[1]):
			continue
		total += 1 # item falls into the predefined class (might not be confident enough!)
	
		genotype = result.genotypeControl (confidence_level = confidence_level) 
		if genotype is None: # enough confidence about the genotype of the healthy cells
			continue
	
		n_test_cases_used += 1

		row, column = 0, 0 
		if genotype == 0.5: 		row = 1
		if genotype == 1.0: 		row = 2
		if result.true_h_vaf == 0.5: 	column = 1
		if result.true_h_vaf == 1.0:	column = 2
		table[row][column] += 1
	
	if n_test_cases_used == 0:
		return None, 0, total, table
	
	return (table[0][0] + table[1][1] + table[2][2]) / float(n_test_cases_used), n_test_cases_used, total, table	


def testCredibleIntervals(results, length_range, confidence_level, cancer_vaf_range):
	n_test_cases_used, total = 0, 0
	in_CI = 0 # number of true cancer VAFs contained in CI

	for result in results:
		if not result.fallsInLengthRange(lower=length_range[0], upper=length_range[1]):
			continue
		if not result.fallsInCancerVAFRange(lower=cancer_vaf_range[0], upper=cancer_vaf_range[1]):
			continue
		total += 1 # item falls into the predefined class (might not be confident enough!)

		#if not result.confidentEnough(confidence_level=confidence_level) or result.global_max_exists == 0:
		#	continue
		n_test_cases_used += 1

		if result.ci[0] <= result.true_c_vaf and result.ci[1] >= result.true_c_vaf:
			in_CI += 1

	if n_test_cases_used == 0:
		return None, 0, total

	return in_CI / float(n_test_cases_used), n_test_cases_used, total	
	
	
def liesInInterval(interval, value):
	if interval[0] != None:
		if value < interval[0]:
			return False
	if interval[1] != None:
		if value > interval[1]:
			return False
	return True
	
def main():

	parser = OptionParser(usage=usage)
	parser.add_option("--applyToAllChromosomes_WARNING_ONLY_IN_THE_END", action="store_true", dest="all_chromosomes", default=True,
						help="Summarize the results for all chromosomes (even and uneven), while normally only the even ones are used. ONLY USE IN THE END! This construction is used to avoid overfitting.")
	parser.add_option("--table", action="store_true", dest="print_tables", default=False,
						help="Prints the 3x3 tables for genotyping and somatic/germline/not present classification.")
	parser.add_option("--truth", action="store_true", dest="use_true_vcf_file", default=False,
						help="Uses the true VCF file to assess the quality of the caller.")
	parser.add_option("--low", action="store_true", dest="low_precision", default=False,
						help="Uses low precision when matching detected indels with the indels in the true VCF file, i.e., the centerpoints must be less or exactly 100 bp away and the length must not differ more than 100 bp. This option overwrites the options '-d' and '-l' (!).")
	parser.add_option("-d", action="store", dest="distance_threshold", default=100, type=int,
				  		help="Distance between centerpoints of deletions threshold. The distance must be smaller than this in order for the deletions to be deemed similar (in addition, the lengths must be similar, see option '-l'). (Default = 100 bp)")
	parser.add_option("-i", action="store", dest="internal_segment_threshold", default=0.0, type=float, 
						help="Number of internal segment alignments from the cancer sample. (Default = 0.0, i.e., no threshold).")
	parser.add_option("-l", action="store", dest="length_threshold", default=100, type=int,
				  		help="Lengths of deletions threshold. The length difference must be smaller than this in order for the deletions to be deemed similar (in addition, the placement of the centerpoints must be similar, see option '-d'). (Default = 100 bp)")
	parser.add_option("-s", action="store", dest="split_read_threshold", default=0.0, type=float, 
						help="Number of split reads from the cancer sample that support the presence of the deletion. (Default = 0.0, i.e., no threshold).")
	parser.add_option("-x", action="store", dest="final_confidence_level", default=0.0, type=float, 
						help="Confidence level used for calling somatic/not somatic. (Default = 0.0, i.e., no threshold).")
	(options, args) = parser.parse_args()
	
	if (len(args)!=1):
		parser.print_help()
		return 1

	if options.low_precision:
		options.length_threshold = 100
		options.distance_threshold = 100

	results, true_c_vafs = obtainResults(open(os.path.abspath(args[0]), 'r'), all_chromosomes = options.all_chromosomes, confidence_level = 0.0, threshold_internal_segment_evidences = options.internal_segment_threshold, threshold_split_read_evidences = options.split_read_threshold)


	list_of_autosomes = range(1,23) 
	#results_new = []
	#for result in results:
	#	if returnIndex(result.chromosome) in list_of_autosomes:
	#		results_new.append(result)
	#results = results_new



	# print genotyping results:
	print('GENOTYPING - CONTROL SAMPLE\n')
	print('length range\tconfidence level\t% correct\t# deletions\ttotal in that range')
	print('------------\t----------------\t---------\t-----------\t-------------------')
	for length_range in LENGTH_RANGES:
		for confidence_level in CONFIDENCE_LEVELS:  
			percentage_correctly_called, n_test_cases_used, total, table = testGenotpyingPerformance(results, length_range, confidence_level)
			
			perc = "--.----"
			if percentage_correctly_called is not None:
				perc = "{:.4f}".format(percentage_correctly_called * 100)
			print(returnIntervalString(length_range), '\t', confidence_level, '\t\t\t', perc, '\t', n_test_cases_used, '\t\t', total)
		print('------------\t----------------\t---------\t-----------\t-------------------')

	if options.print_tables:
		for length_range in LENGTH_RANGES:
			for confidence_level in CONFIDENCE_LEVELS:  
				percentage_correctly_called, n_test_cases_used, total, table = testGenotpyingPerformance(results, length_range, confidence_level)
				print('')
				title = "GENOTYPING " +  returnTitleTable(length_range, confidence_level, CANCER_VAF_RANGE, None)	
				printTable3x3(table, title = title, labels = ["0.0", "0.5", "1.0"])
				print('')
				

	# print main test results 
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
	
	if options.print_tables:
		for length_range in LENGTH_RANGES:
			for confidence_level in CONFIDENCE_LEVELS:  
				percentage_correctly_called, n_test_cases_used, total, table = testCallingPerformance(results, length_range, confidence_level, CANCER_VAF_RANGE)
				print('')
				title = returnTitleTable(length_range, confidence_level, CANCER_VAF_RANGE, None)	
				printTable3x3(table, title = title, labels = ["som.", "germl", "n.prs"])
				print('')

	print('\nSOMATIC/NOT SOMATIC - germline and not present are now considered to be one class\n')
	print('length range\tconfidence level\trecall\t\tprecision\t# deletions\ttotal in that range')
	print('------------\t----------------\t------\t\t---------\t-----------\t-------------------')
	for length_range in LENGTH_RANGES:
		for confidence_level in CONFIDENCE_LEVELS:  
			recall, precision, n_test_cases_used, total = testSomaticMutationCallingPerformance(results, length_range, confidence_level, CANCER_VAF_RANGE)
			recall_string, precision_string = '-.----', '-.----'				
			if recall is not None:
				recall_string = "{:.4f}".format(recall)
			if precision is not None:
				precision_string = "{:.4f}".format(precision)
			print(returnIntervalString(length_range), '\t', confidence_level, '\t\t\t', recall_string, '\t', precision_string, '\t', n_test_cases_used, '\t\t', total)
		print('------------\t----------------\t------\t\t---------\t-----------\t-------------------')


	print('\nPRESENT/NOT PRESENT - can the caller distinguish between an indel being present or not?\n')
	print('length range\tconfidence level\trecall\t\tprecision\t# deletions\ttotal in that range')
	print('------------\t----------------\t------\t\t---------\t-----------\t-------------------')
	for length_range in LENGTH_RANGES:
		for confidence_level in CONFIDENCE_LEVELS:  
			recall, precision, n_test_cases_used, total = testPresenceDetectionPerformance(results, length_range, confidence_level, CANCER_VAF_RANGE)
			recall_string, precision_string = '-.----', '-.----'				
			if recall is not None:
				recall_string = "{:.4f}".format(recall)
			if precision is not None:
				precision_string = "{:.4f}".format(precision)
			print(returnIntervalString(length_range), '\t', confidence_level, '\t\t\t', recall_string, '\t', precision_string, '\t', n_test_cases_used, '\t\t', total)
		print('------------\t----------------\t------\t\t---------\t-----------\t-------------------')

	print('\nTEST CREDIBLE INTERVALS\n')
	print('confidence interval\tlength range\tconfidence level\t% true values in CI\t# deletions\ttotal in that range')
	print('-------------------\t------------\t----------------\t-------------------\t-----------\t-------------------')
	for length_range in LENGTH_RANGES:
		for confidence_level in CONFIDENCE_LEVELS:  
			in_CI, n_test_cases_used, total = testCredibleIntervals(results, length_range, confidence_level, CANCER_VAF_RANGE)
			in_CI_string = "--.----"
			if in_CI is not None:	
				in_CI_string = "{:.4f}".format(in_CI * 100)
			print('95%', '\t\t\t', returnIntervalString(length_range), '\t', confidence_level, '\t\t\t', in_CI_string, '\t\t', n_test_cases_used, '\t\t', total)
		print('-------------------\t------------\t----------------\t-------------------\t-----------\t-------------------')

	print('\nMEAN SQUARED ERROR (CANCER VAF)\n')
	print('length range\tconfidence level\tmean squared error')
	print('------------\t----------------\t------------------')
	for length_range in LENGTH_RANGES:
		for confidence_level in CONFIDENCE_LEVELS:
			mean_squared_error = 0.0
			n_obs = 0.0 
			for result in results:
				if not result.fallsInLengthRange(length_range[0], length_range[1]):
					continue	
				if not result.confidentEnoughAboutCalling(confidence_level = confidence_level):
					continue
				if result.mle_c_vaf == None:
					continue
				mean_squared_error += (result.true_c_vaf - result.mle_c_vaf)**2.0
				n_obs += 1.0
			mse_str = '-.----'
			if n_obs != 0.0:
				mse_str = "{:.4f}".format(mean_squared_error / n_obs)	
			print(returnIntervalString(length_range), '\t', confidence_level, '\t\t\t', mse_str)
		print('------------\t----------------\t------------------')	


	if options.use_true_vcf_file: # use the true VCF file to assess the quality of the caller at hand
		print('\nReading in the true deletions from the true VCF file...')
		obtainListOfTrueDeletions(10) # read in all the true deletions from the true VCF file
		print('Done reading in the true deletions from the true VCF file...\n')
		
		print('cancer vaf range\tlength range\t# deletions\t# somatic\t# called somatic\t# true positives\tRECALL\tPRECISION')
		print('----------------\t------------\t-----------\t---------\t----------------\t----------------\t------\t---------')

		for cancer_vaf_range in CANCER_VAF_RANGES:
			# initialize data structures for storing the results
			n_true_deletions 	= [0 for i in range(len(LENGTH_RANGES))] # number of true deletions per length range
			n_true_somatic 		= [0 for i in range(len(LENGTH_RANGES))] # number of true somatic deletions per length range
			n_called_somatic 	= [0 for i in range(len(LENGTH_RANGES))] # number of deletions that we detected and called somatic
			n_true_positives 	= [0 for i in range(len(LENGTH_RANGES))] # number of deletions that we called somatic and were indeed somatic 
		
			# determine the number of calls we made and how many were correct
			for result in results:
				call 	= result.call(confidence_level = options.final_confidence_level)  
				truth 	= result.obtainTruth()

				if call == 'None':
					continue
	
				if not result.fallsInCancerVAFRange(lower = cancer_vaf_range[0], upper = cancer_vaf_range[1]):
					continue

				for i in range(len(LENGTH_RANGES)):
					length_range = LENGTH_RANGES[i] 	
					if result.fallsInLengthRange(lower = length_range[0], upper = length_range[1]):
						if call == 'somatic':
							n_called_somatic[i] += 1
						if call == 'somatic' and truth == 'somatic':
							n_true_positives[i] += 1

			for autosome in list_of_autosomes: 
				for vcf_record in list_true_deletions[autosome]:
					true_deletion = Deletion(vcf_record)
					is_somatic = somatic(80, vcf_record)
					if is_somatic:
						true_h_vaf, true_c_vaf = returnTrueVAFs (80, vcf_record)
						if not liesInInterval(cancer_vaf_range, true_c_vaf):
							continue

					for i in range(len(LENGTH_RANGES)):
						length_range = LENGTH_RANGES[i]
						if true_deletion.fallsInLengthRange(length_range[0], length_range[1]):
							n_true_deletions[i] += 1
							if is_somatic:
								n_true_somatic[i] += 1

			for i in range(len(LENGTH_RANGES)):
				print(returnIntervalString(cancer_vaf_range), '\t\t', returnIntervalString(LENGTH_RANGES[i]), '\t', n_true_deletions[i], '\t\t', n_true_somatic[i], '\t\t', n_called_somatic[i], '\t\t\t', n_true_positives[i], '\t\t\t', end = '')
				recall_str = '-.----'
				precision_str = '-.----'
				if n_true_somatic[i] != 0:
					recall_str = "{:.4f}".format(float(n_true_positives[i]) / float(n_true_somatic[i]))
				if n_called_somatic[i] != 0:
					precision_str = "{:.4f}".format(float(n_true_positives[i]) / float(n_called_somatic[i]))				
				print(recall_str, '\t', precision_str)
			print('----------------\t------------\t-----------\t---------\t----------------\t----------------\t------\t---------')
			#for i in range(len(LENGTH_RANGES)):
			#	print("LENGTH RANGE: ", returnIntervalString(LENGTH_RANGES[i]))	
			#	print('total number of true deletions in the VCF file: ', n_true_deletions[i])
			#	print('total number of true somatic deletions in the VCF file: ', n_true_somatic[i])
			#	print('total number of somatic predictions made: ', n_called_somatic[i])
			#	print('true positives: ', n_true_positives[i])
			#	print('recall: ', float(n_true_positives[i]) / float(n_true_somatic[i]))
			#	print('precision: ', float(n_true_positives[i]) / float(n_called_somatic[i]))		
			#	print('')


if __name__ == '__main__':
	sys.exit(main())
