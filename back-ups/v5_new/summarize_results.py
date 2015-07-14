#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
from scipy.stats import norm 
from scipy.stats import chi2 
import os
import sys
import math

from ResultsInstance import *

__author__ = "Louis Dijkstra"

usage = """%prog <results-filename>

	<results-filename> 	File with the somatic mutation results. 
	
NOTE: 	At the moment we only process cases where the cancer VAF is assumed
		continuous. In addition, only autosomes are taken into account.
"""

LENGTH_RANGES = [[None, None], [10,29], [30,49], [50,59], [60,69], [70,79],[80,89],[90,99]]
CANCER_VAF_RANGES = [[None, None]]#[[None, None], [0.0, 0.20], [.20, .40], [.40, .60], [.60, .80], [.80, 1.0]]
CONFIDENCE_LEVELS = [0.0, .90, .95, .99]
CONFIDENCE_INTERVALS = [.90, .95]
CANCER_VAF_RANGE = [None, None]

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

def testCallingPerformance(results, length_range, confidence_level, cancer_vaf_range, confidence_interval, p_value_threshold, max_width_ci, use_ci = False):
	"""Returns the performance for calling an indel somatic/germline or not present."""
	n_test_cases_used, total = 0, 0
	table = [[0,0,0],[0,0,0],[0,0,0]]
	for result in results:
		if not result.fallsInLengthRange(lower=length_range[0], upper=length_range[1]):
			continue
		if not result.fallsInCancerVAFRange(lower=cancer_vaf_range[0], upper=cancer_vaf_range[1]):
			continue
		total += 1 # item falls into the predefined class (might not be confident enough!)
		
		call = result.call(p_value_threshold = p_value_threshold, confidence_level = confidence_level, confidence_interval = confidence_interval, width_conf_interval = max_width_ci, use_ci = use_ci)
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

def testSomaticMutationCallingPerformance(results, length_range, confidence_level, cancer_vaf_range, confidence_interval, p_value_threshold, max_width_ci, use_ci = False):
	n_test_cases_used, total = 0, 0
	table = [[0,0],[0,0],[0,0]]
	for result in results:
		if not result.fallsInLengthRange(lower=length_range[0], upper=length_range[1]):
			continue
		if not result.fallsInCancerVAFRange(lower=cancer_vaf_range[0], upper=cancer_vaf_range[1]):
			continue
		total += 1 # item falls into the predefined class (might not be confident enough!)

		call = result.call(p_value_threshold = p_value_threshold, confidence_level = confidence_level, confidence_interval = confidence_interval, width_conf_interval = max_width_ci, use_ci = use_ci)
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


def testPresenceDetectionPerformance(results, length_range, confidence_level, cancer_vaf_range, confidence_interval, p_value_threshold, max_width_ci, use_ci = False):
	n_test_cases_used, total = 0, 0
	table = [[0,0],[0,0],[0,0]]
	for result in results:
		if not result.fallsInLengthRange(lower=length_range[0], upper=length_range[1]):
			continue
		if not result.fallsInCancerVAFRange(lower=cancer_vaf_range[0], upper=cancer_vaf_range[1]):
			continue
		total += 1 # item falls into the predefined class (might not be confident enough!)

		call = result.call(p_value_threshold = p_value_threshold, confidence_level = confidence_level, confidence_interval = confidence_interval, width_conf_interval = max_width_ci, use_ci = use_ci)
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


def testConfidenceIntervals(results, length_range, confidence_level, cancer_vaf_range, confidence_interval):
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

		ci = None
		if confidence_interval == .95:
			ci = result.ci_95_0
		if confidence_interval == .90:
			ci = result.ci_90_0

		if result.maxlogl1 > result.maxlogl0 and result.maxlogl1 > result.maxlogl2:
			if confidence_interval == .95:
				ci = result.ci_95_1
			if confidence_interval == .90:
				ci = result.ci_90_1
		if result.maxlogl2 > result.maxlogl0 and result.maxlogl2 > result.maxlogl1:
			if confidence_interval == .95:
				ci = result.ci_95_2
			if confidence_interval == .90:
				ci = result.ci_90_2

		n_test_cases_used += 1
		
		if ci[0] <= result.true_c_vaf and ci[1] >= result.true_c_vaf:
			in_CI += 1

	if n_test_cases_used == 0:
		return None, 0, total

	return in_CI / float(n_test_cases_used), n_test_cases_used, total	
	
		
def main():

	parser = OptionParser(usage=usage)
	parser.add_option("--applyToAllChromosomes_WARNING_ONLY_IN_THE_END", action="store_true", dest="all_chromosomes", default=False,
						help="Summarize the results for all chromosomes (even and uneven), while normally only the even ones are used. ONLY USE IN THE END! This construction is used to avoid overfitting.")
	parser.add_option("--bonf", action="store_true", dest="bonferroni", default=False,
						help="A Bonferroni correction is applied to the p-value threshold.")
	parser.add_option("--table", action="store_true", dest="print_tables", default=False,
						help="Prints the 3x3 tables for genotyping and somatic/germline/not present classification.")
	parser.add_option("--use-ci", action="store_true", dest="use_ci", default=False,
						help="The 95% confidence interval is used to make the call, i.e., when 0 is in the 95% CI, we call it 'not present', otherwise 'somatic'.")
	parser.add_option("-i", action="store", dest="internal_segment_threshold", default=0.0, type=float, 
						help="Number of internal segment alignments from the cancer sample. (Default = 0.0, i.e., no threshold).")
	parser.add_option("-p", action="store", dest="p_value_threshold", default=0.05, type=float, 
						help="P-value threshold used for deciding whether an indel is somatic/not present. (Default = 0.05).")
	parser.add_option("-s", action="store", dest="split_read_threshold", default=0.0, type=float, 
						help="Number of split reads from the cancer sample that support the presence of the deletion. (Default = 0.0, i.e., no threshold).")
	parser.add_option("-w", action="store", dest="max_width_ci", default=1.0, type=float,
						help="Threshold on the width of the confidence interval. (Default = 1.0; all results are taken into account)")
	(options, args) = parser.parse_args()
	
	if (len(args)!=1):
		parser.print_help()
		return 1

	results, true_c_vafs = obtainResults(open(os.path.abspath(args[0]), 'r'), all_chromosomes = options.all_chromosomes, confidence_level = 0.0, threshold_internal_segment_evidences = options.internal_segment_threshold, threshold_split_read_evidences = options.split_read_threshold)

	p_value_threshold = options.p_value_threshold 
	if options.bonferroni:
		p_value_threshold /= float(len(results))

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
	print('length range\tconfidence level\tconfidence interval\t% correct\t# deletions\ttotal in that range')
	print('------------\t----------------\t-------------------\t---------\t-----------\t-------------------')
	for length_range in LENGTH_RANGES:
		for confidence_level in CONFIDENCE_LEVELS:  
			for confidence_interval in CONFIDENCE_INTERVALS:
				percentage_correctly_called, n_test_cases_used, total, table = testCallingPerformance(results, length_range, confidence_level, CANCER_VAF_RANGE, confidence_interval, p_value_threshold, options.max_width_ci, use_ci = options.use_ci)	
				perc = "--.----"
				if percentage_correctly_called is not None:				
					perc = "{:.4f}".format(percentage_correctly_called * 100)
				print(returnIntervalString(length_range), '\t', confidence_level, '\t\t\t', confidence_interval,'\t\t\t', perc, '\t', n_test_cases_used, '\t\t', total)
		print('------------\t----------------\t-------------------\t---------\t-----------\t-------------------')
	
	if options.print_tables:
		for length_range in LENGTH_RANGES:
			for confidence_level in CONFIDENCE_LEVELS:  
				for confidence_interval in CONFIDENCE_INTERVALS:
					percentage_correctly_called, n_test_cases_used, total, table = testCallingPerformance(results, length_range, confidence_level, CANCER_VAF_RANGE, confidence_interval, p_value_threshold, options.max_width_ci, use_ci = options.use_ci)
					print('')
					title = returnTitleTable(length_range, confidence_level, CANCER_VAF_RANGE, None)	
					printTable3x3(table, title = title, labels = ["som.", "germl", "n.prs"])
					print('')

	print('\nSOMATIC/NOT SOMATIC - germline and not present are now considered to be one class\n')
	print('length range\tconfidence level\tconfidence interval\trecall\t\tprecision\t# deletions\ttotal in that range')
	print('------------\t----------------\t-------------------\t------\t\t---------\t-----------\t-------------------')
	for length_range in LENGTH_RANGES:
		for confidence_level in CONFIDENCE_LEVELS:  
			for confidence_interval in CONFIDENCE_INTERVALS:
				recall, precision, n_test_cases_used, total = testSomaticMutationCallingPerformance(results, length_range, confidence_level, CANCER_VAF_RANGE, confidence_interval, p_value_threshold, options.max_width_ci, use_ci = options.use_ci)
				recall_string, precision_string = '-.----', '-.----'				
				if recall is not None:
					recall_string = "{:.4f}".format(recall)
				if precision is not None:
					precision_string = "{:.4f}".format(precision)
				print(returnIntervalString(length_range), '\t', confidence_level, '\t\t\t', confidence_interval,'\t\t\t', recall_string, '\t', precision_string, '\t', n_test_cases_used, '\t\t', total)
		print('------------\t----------------\t-------------------\t------\t\t---------\t-----------\t-------------------')


	print('\nPRESNENT/NOT PRESENT - can the caller distinguish between an indel being present or not?\n')
	print('length range\tconfidence level\tconfidence interval\trecall\t\tprecision\t# deletions\ttotal in that range')
	print('------------\t----------------\t-------------------\t------\t\t---------\t-----------\t-------------------')
	for length_range in LENGTH_RANGES:
		for confidence_level in CONFIDENCE_LEVELS:  
			for confidence_interval in CONFIDENCE_INTERVALS:
				recall, precision, n_test_cases_used, total = testPresenceDetectionPerformance(results, length_range, confidence_level, CANCER_VAF_RANGE, confidence_interval, p_value_threshold, options.max_width_ci, use_ci = options.use_ci)
				recall_string, precision_string = '-.----', '-.----'				
				if recall is not None:
					recall_string = "{:.4f}".format(recall)
				if precision is not None:
					precision_string = "{:.4f}".format(precision)
				print(returnIntervalString(length_range), '\t', confidence_level, '\t\t\t', confidence_interval,'\t\t\t', recall_string, '\t', precision_string, '\t', n_test_cases_used, '\t\t', total)
		print('------------\t----------------\t-------------------\t------\t\t---------\t-----------\t-------------------')

	print('\nTEST CONFIDENCE INTERVALS\n')
	print('confidence interval\tlength range\tconfidence level\t% true values in CI\t# deletions\ttotal in that range')
	print('-------------------\t------------\t----------------\t-------------------\t-----------\t-------------------')
	for confidence_interval in CONFIDENCE_INTERVALS:
		for length_range in LENGTH_RANGES:
			for confidence_level in CONFIDENCE_LEVELS:  
				in_CI, n_test_cases_used, total = testConfidenceIntervals(results, length_range, confidence_level, CANCER_VAF_RANGE, confidence_interval)
				in_CI_string = "{:.4f}".format(in_CI * 100)
				print(confidence_interval, '\t\t\t', returnIntervalString(length_range), '\t', confidence_level, '\t\t\t', in_CI_string, '\t\t', n_test_cases_used, '\t\t', total)
		print('-------------------\t------------\t----------------\t-------------------\t-----------\t-------------------')



				
	
if __name__ == '__main__':
	sys.exit(main())
