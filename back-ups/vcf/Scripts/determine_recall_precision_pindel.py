#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
import vcf
import os
import sys
import math

from Alignments import * 
from Indel import *

__author__ = "Louis Dijkstra"

usage = """%prog options <delly-vcf>

Determines the recall and precision (as defined in the MATE-CLEVER paper) for the deletions in the raw PINDEL VCF. 

	<pindel-vcf> 	The PINDEL VCF file as provided by Tobias.
"""

LENGTH_RANGES 		= [[None, None], [10,29], [30,49], [50,69], [70,99], [100, 249], [250, None]]
CANCER_VAF_RANGES 	= [[None, None], [0.0, 0.20], [.20, .40], [.40, .60], [.60, .80], [.80, 1.0], [.20, None], [.40, None], [.50, None], [.60, None]]

VCF_TRUTH_FILENAME = "/data1/structvar/vaf-experiments/vcf/cancerclones.vcf" # VCF file that contains the ground truth
VCF_TRUTH_READER = vcf.Reader(open(VCF_TRUTH_FILENAME)) # used for reading the VCF file 

list_true_deletions = [[] for i in range(24)] # 24 empty lists for the 22 autosomes and the two sex-chromosomes X and Y   

class DIndel: 
	"""Class for representing an indel (deletion/insertion)."""
	def __init__(self, vcf_record):
		self.vcf_record = vcf_record 
		self.chromosome = str(vcf_record.CHROM)
		self.length 	= returnLength(vcf_record)	

	def fallsInLengthRange (self, lower = None, upper = None):
		"""Checks whether an indel is in a length range. In case of 'None', that side of the interval is considered unbounded, e.g., [250, None] are all indels of length larger than 250"""
		if lower != None:
			if self.length < lower:
				return False
		if upper != None:
			if self.length > upper:
				return False
		return True
	
class DDeletion(DIndel):
	"""Class for representing a deletion. """
	def __init__(self, vcf_record):
		DIndel.__init__(self, vcf_record)
		# Interval which is deleted [start, end]:
		self.start = self.vcf_record.POS + 1
		self.end   = self.start + self.length - 1
		# centerpoints are the points in the middle of the deletion; when length is even, there are two, otherwise one
		if self.length % 2 == 0: 
			self.centerpoints = [self.start + self.length / 2 - 1, self.start + self.length / 2]
		else: 
			self.centerpoints = [self.start + self.length / 2]

	def print(self):
		"""Prints most relevant data on the deletion to standard output"""
		print(self.chromosome, '\t', self.vcf_record.POS, '\t', self.length, '\t', end = "")

	def similarTo(self, other_deletion, difference_length_threshold = 100, difference_centerpoint_threshold = 100):
		"""Returns True when the other_deletion is considered similar. Otherwise False."""
		if returnMinimumDifference(self.centerpoints, other_deletion.centerpoints) > difference_centerpoint_threshold:
			return False
		if abs(self.length - other_deletion.length) > difference_length_threshold:
			return False
		return True


def strip(string):
	string = string.lstrip()
	return string.rstrip()

class DellyDeletion(DDeletion):
	"""Used for representing a Deletion from the Delly VCF."""
	def __init__(self, vcf_record):
		DDeletion.__init__(self, vcf_record)
		self.is_somatic = self.isSomatic()		

	def isSomatic(self):
		control, cancer = None, None
		for call in self.vcf_record.samples:
			if call.sample == 'NORMAL':
				if call.gt_nums != None: 
					control = strip(call.gt_nums)
			if call.sample == 'CANCER':
				if call.gt_nums != None: 				
					cancer = strip(call.gt_nums)	
		if control == '0/0':
			if cancer != '0/0' and cancer != './.':
				return True	
		return False		
	
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

def returnIntervalString (interval):
	if interval[0] is None and interval[1] is not None:
		return "at most " + str(interval[1])
	elif interval[0] is not None and interval[1] is None:
		return "at least " + str(interval[0])
	elif interval[0] is None and interval[1] is None:
		return "unbounded"
	else:
		return "[" + str(interval[0]) + ','  + str(interval[1]) + "]"

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


def liesInInterval(interval, value):
	if interval[0] != None:
		if value < interval[0]:
			return False
	if interval[1] != None:
		if value > interval[1]:
			return False
	return True

def isPINDELDeletion(pindel_vcf_record):
	if 'SVTYPE' in pindel_vcf_record.INFO:
		if pindel_vcf_record.INFO['SVTYPE'] == 'DEL':
			return True
	return False

def returnLength(pindel_vcf_record):
	return abs(int(pindel_vcf_record.INFO['SVLEN']))

def main():

	parser = OptionParser(usage=usage)
	parser.add_option("--low", action="store_true", dest="low_precision", default=False,
						help="Uses low precision when matching detected indels with the indels in the true VCF file, i.e., the centerpoints must be less or exactly 100 bp away and the length must not differ more than 100 bp. This option overwrites the options '-d' and '-l' (!).")
	parser.add_option("-d", action="store", dest="distance_threshold", default=100, type=int,
				  		help="Distance between centerpoints of deletions threshold. The distance must be smaller than this in order for the deletions to be deemed similar (in addition, the lengths must be similar, see option '-l'). (Default = 100 bp)")
	parser.add_option("-l", action="store", dest="length_threshold", default=100, type=int,
				  		help="Lengths of deletions threshold. The length difference must be smaller than this in order for the deletions to be deemed similar (in addition, the placement of the centerpoints must be similar, see option '-d'). (Default = 100 bp)")
	(options, args) = parser.parse_args()


	list_of_autosomes = range(1,2) # TODO TODO change! 

	if (len(args)!=1):
		parser.print_help()
		return 1

	if options.low_precision:
		options.length_threshold = 100
		options.distance_threshold = 100
	
	delly_vcf_reader = vcf.Reader(open(args[0]))
	somatic_calls = [[] for i in range(24)] # all somatic deletions from Delly are stored per chromosome

	print('Reading in the PINDEL VCF...')
	n_called_somatic = [0 for i in range(len(LENGTH_RANGES))] # number of deletions that we detected and called somatic
	for vcf_record in delly_vcf_reader:  
		if not isPINDELDeletion(vcf_record):
			continue
		if returnLength(vcf_record) < 10:
			continue
		deletion = DellyDeletion(vcf_record)
		if not returnIndex(deletion.chromosome) in list_of_autosomes:
			continue
		if deletion.is_somatic: 
			somatic_calls[returnIndex(deletion.chromosome)].append(deletion)
			for i in range(len(LENGTH_RANGES)):
				length_range = LENGTH_RANGES[i]
				if deletion.fallsInLengthRange(lower=length_range[0], upper=length_range[1]):
					n_called_somatic[i] += 1
	print('Done reading in the PINDEL VCF...')
	
	print('Reading true VCF file @ /data1/structvar/vaf-experiments/vcf/cancerclones.vcf')
	obtainListOfTrueDeletions(10) # read in all the true deletions from the true VCF file	
	print('Done reading in the true deletions from the true VCF file...\n')
	
	n_true_deletions 	= [0 for i in range(len(LENGTH_RANGES))] # number of true deletions per length range
	n_true_somatic 		= [0 for i in range(len(LENGTH_RANGES))] # number of true somatic deletions per length range
	n_true_positives 	= [0 for i in range(len(LENGTH_RANGES))] # number of deletions that we called somatic and were indeed somatic 
	for autosome in list_of_autosomes:
		for true_vcf_record in list_true_deletions[autosome]:
			true_deletion = Deletion(true_vcf_record) 
			is_somatic = somatic(80, true_vcf_record)
			called_correctly = False
			if is_somatic: # determine whether deletion was called correctly or not 
				for delly_deletion in somatic_calls[autosome]:
					if returnMinimumDifference(delly_deletion.centerpoints, true_deletion.centerpoints) > options.distance_threshold:
						continue
					if abs(delly_deletion.length - true_deletion.length) > options.length_threshold:
						continue
					called_correctly = True

			#print('TRUE DELETION')
			#true_deletion.print()
			#print('\tlength', true_deletion.length)
			#print('')
			#print(is_somatic, '\t', called_correctly)
			for i in range(len(LENGTH_RANGES)):
				length_range = LENGTH_RANGES[i]
				if true_deletion.fallsInLengthRange(length_range[0], length_range[1]):
					#print(returnIntervalString(length_range), ' deletion falls in here!')
					n_true_deletions[i] += 1
					if is_somatic: 
						n_true_somatic[i] += 1
					if called_correctly:
						n_true_positives[i] += 1

	print('cancer vaf range\tlength range\t# deletions\t# somatic\t# called somatic\t# true positives\tRECALL\tPRECISION')
	print('----------------\t------------\t-----------\t---------\t----------------\t----------------\t------\t---------')
	for i in range(len(LENGTH_RANGES)):
		print(returnIntervalString([None, None]), '\t\t', returnIntervalString(LENGTH_RANGES[i]), '\t', n_true_deletions[i], '\t\t', n_true_somatic[i], '\t\t', n_called_somatic[i], '\t\t\t', n_true_positives[i], '\t\t\t', end = '')
		recall_str = '-.-----'
		precision_str = '-.-----'
		if n_true_somatic[i] != 0:
			recall_str = "{:.4f}".format(float(n_true_positives[i]) / float(n_true_somatic[i]))
		if n_called_somatic[i] != 0:
			precision_str = "{:.4f}".format(float(n_true_positives[i]) / float(n_called_somatic[i]))				
		print(recall_str, '\t', precision_str)
	print('----------------\t------------\t-----------\t---------\t----------------\t----------------\t------\t---------')

	print('')

	print('cancer vaf range\tlength range\t# deletions\t# somatic\t# called somatic\t# true positives\tRECALL\tPRECISION')
	print('----------------\t------------\t-----------\t---------\t----------------\t----------------\t------\t---------')
	for cancer_vaf_range in CANCER_VAF_RANGES:
		n_true_deletions 	= [0 for i in range(len(LENGTH_RANGES))] # number of true deletions per length range
		n_true_somatic 		= [0 for i in range(len(LENGTH_RANGES))] # number of true somatic deletions per length range
		n_true_positives 	= [0 for i in range(len(LENGTH_RANGES))] # number of deletions that we called somatic and were indeed somatic 
		for autosome in list_of_autosomes:
			for true_vcf_record in list_true_deletions[autosome]:
				true_deletion = Deletion(true_vcf_record) 
				true_h_vaf, true_c_vaf = returnTrueVAFs (80, true_vcf_record)
				if not liesInInterval(cancer_vaf_range, true_c_vaf):
					continue
				is_somatic = somatic(80, true_vcf_record)
				called_correctly = False
				if is_somatic: # determine whether deletion was called correctly or not 
					for delly_deletion in somatic_calls[autosome]:
						if returnMinimumDifference(delly_deletion.centerpoints, true_deletion.centerpoints) > options.distance_threshold:
							continue
						if abs(delly_deletion.length - true_deletion.length) > options.length_threshold:
							continue
						called_correctly = True
				for i in range(len(LENGTH_RANGES)):
					length_range = LENGTH_RANGES[i]
					if true_deletion.fallsInLengthRange(length_range[0], length_range[1]):
						n_true_deletions[i] += 1
						if is_somatic: 
							n_true_somatic[i] += 1
						if called_correctly:
							n_true_positives[i] += 1


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
"""




	true_somatic = [[] for i in range(24)] # stores booleans; true when true deletion is somatic, false otherwise.
	called_correctly = [[] for i in range(24)] # stores for every true deletion a boolean whether it is called correctly
	for autosome in range(1,23):
		for vcf_record in list_true_deletions[autosome]:
			som, cc = False, False # somatic/called correctly
			true_deletion = Deletion(vcf_record) 
			if somatic(80, vcf_record):
				som = True
				for delly_deletion in somatic_calls[autosome]: # walk through all deletions delly thought were somatic
					if returnMinimumDifference(delly_deletion.centerpoints, true_deletion.centerpoints) > options.distance_threshold:
						continue
					if abs(delly_deletion.length - true_deletion.length) > options.length_threshold:
						continue
					cc = True	
			true_somatic[autosome].append(som)
			called_correctly[autosome].append(cc)
				
	print('cancer vaf range\tlength range\t# deletions\t# somatic\t# called somatic\t# true positives\tRECALL\tPRECISION')
	print('----------------\t------------\t-----------\t---------\t----------------\t----------------\t------\t---------')

	for cancer_vaf_range in CANCER_VAF_RANGES:
		# initialize data structures for storing the results
		n_true_deletions 	= [0 for i in range(len(LENGTH_RANGES))] # number of true deletions per length range
		n_true_somatic 		= [0 for i in range(len(LENGTH_RANGES))] # number of true somatic deletions per length range
		n_called_somatic 	= [0 for i in range(len(LENGTH_RANGES))] # number of deletions that we detected and called somatic
		n_true_positives 	= [0 for i in range(len(LENGTH_RANGES))] # number of deletions that we called somatic and were indeed somatic 
		
		for autosome in range(1,23): 
			for i in range(len(list_true_deletions[autosome])):
				vcf_record = list_true_deletions[autosome][i]
				true_deletion = Deletion(vcf_record)
				true_h_vaf, true_c_vaf = returnTrueVAFs (80, vcf_record)
				som = true_somatic[autosome][i]
				cc = called_correctly[autosome][i]
				if som and not liesInInterval(cancer_vaf_range, true_c_vaf):
					continue
				for j in range(len(LENGTH_RANGES)):
					length_range = LENGTH_RANGES[j]
					if true_deletion.fallsInLengthRange(lower = length_range[0], upper = length_range[1]):
						n_true_deletions[j] += 1
						if som: n_true_somatic[j] += 1
						if cc: n_true_positives[j] += 1

		# determine called somatic
		for autosome in range(1,23):
			for delly_deletion in somatic_calls[autosome]:
				for j in range(len(LENGTH_RANGES)):
					length_range = LENGTH_RANGES[j]
					if delly_deletion.fallsInLengthRange(lower = length_range[0], upper = length_range[1]):	
						n_called_somatic[j] += 1
				

		for i in range(len(LENGTH_RANGES)):
			print(returnIntervalString(cancer_vaf_range), '\t\t', returnIntervalString(LENGTH_RANGES[i]), '\t', n_true_deletions[i], '\t\t', n_true_somatic[i], '\t\t', n_called_somatic[i], '\t\t\t', n_true_positives[i], '\t\t', end = '')
			recall_str = '--.-----'
			precision_str = '--.-----'
			if n_true_somatic[i] != 0:
				recall_str = "{:.4f}".format(float(n_true_positives[i]) / float(n_true_somatic[i]))
			if n_called_somatic[i] != 0:
				precision_str = "{:.4f}".format(float(n_true_positives[i]) / float(n_called_somatic[i]))				
			print(recall_str, '\t', precision_str)
		print('----------------\t------------\t-----------\t---------\t----------------\t----------------\t------\t---------')
"""
if __name__ == '__main__':
	sys.exit(main())




