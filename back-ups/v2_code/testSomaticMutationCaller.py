#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
from collections import defaultdict
from scipy.stats import norm 
from scipy.optimize import fminbound
from scipy.optimize import newton
import os
import sys
import vcf
import pysam
import math


# import from SomaticMutationCaller.py
from SomaticMutationCaller import * 


# import from VAFEstimator.py
from VAFEstimator import *



__author__ = "Louis Dijkstra"

usage = """%prog [options] <mu> <sigma> <vcf-filename> <bam-healthy> <bam-cancer> 

	Tests the somatic mutation caller

	<mu>		the mean of the null normal 
				distribution for the internal 
				segment length when not 
				affected by an indel. 
	<sigma>		the standard deviation of the 
				null normal distribution for 
				the internal segment length 
				when not affected by an indel
	<vcf-file> 	tabix-indexed VCF file containing the 
				indels which allele frequency 
				needs to be estimated. 
				NOTE: File needs to be sorted 
				and indexed. 
	<bam-healthy>	BAM file containing the 
				alignments of the disease sample. 
				NOTE: File needs to be sorted 
				and indexed. 
	<bam-cancer>	BAM file containing the 
				alignments of the disease sample. 
				NOTE: File needs to be sorted 
				and indexed. 


NOTE: the current implementation
		1) 	ignores indels in the vcf file 
			for which more than one 
			alternative (ALT) are provided. 	
		2)	only considers primary 
			alignments; secondary/ternary etc. 
			alignments are neglected. 
"""



class SMTestCase:
	"""
		Class that contains the relevant information for one indel.
		Used for testing the somatic mutation caller.
	"""

	def __init__(self, is_deletion, exists, h_true_vaf, c_true_vaf, true_sm, call, h_mle, c_mle, ci, unique, indel, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments):
		self.is_deletion 		= is_deletion
		self.exists			= exists 
		self.h_true_vaf 		= h_true_vaf
		self.c_true_vaf			= c_true_vaf
		self.true_sm			= true_sm
		self.call 			= call 		
		self.h_mle 			= h_mle
		self.c_mle			= c_mle
		self.ci				= ci 
		self.unique			= unique
		self.indel			= indel
		self.position			= indel.vcf_record.POS
		self.length			= indel.length
		self.n_h_paired_end_alignments 	= len(h_paired_end_alignments)
		self.n_h_overlapping_alignments = len(h_overlapping_alignments)
		self.n_c_paired_end_alignments 	= len(c_paired_end_alignments)
		self.n_c_overlapping_alignments = len(c_overlapping_alignments)

		if self.true_sm == self.call:
			self.correct = True
		else:
			self.correct = False 
		
		self.Xh, self.Yh, self.Xc, self.Yc = [], [], [], []
	
		for a in h_paired_end_alignments: 
			self.Xh.append(a.value)

		for a in h_overlapping_alignments:
			self.Yh.append(a.value)

		for a in c_paired_end_alignments: 
			self.Xc.append(a.value)

		for a in c_overlapping_alignments:
			self.Yc.append(a.value)

	def print(self, min_length = 0):
		"""Prints relevant information to the command line"""
		if self.length >= min_length: 
			print(self.length, '\t|\t', self.exists, '\t', self.correct, '\t|\t', self.true_sm, '\t',self.call, '\t|\t', self.h_true_vaf, self.c_true_vaf, '\t|\t', self.h_mle, self.c_mle, self.ci)

	def isOfLength(self, a, b):
		"""Determines whether this test case falls into a certain length class"""
		if self.length >= a and self.length <= b:
			return True
		return False




class SMTestDataCollector: 
	"""
		Collects/stores data used for testing
	"""
	def __init__(self, truth = False):
		self.testcases = [] 
		self.length_range = [[10, 29],[30,49],[50,99],[100,249]]
		self.min_length = 10 
		self.truth = truth 
		self.true_vcf_reader = vcf.Reader(open("/data1/structvar/vaf-experiments/vcf/cancerclones.vcf.gz")) 


		if truth:
			print("Only the VCF file with the true indels will be used for analysis.") 
			# Get the true cancer VCF file
			

	def determineVAF (self, gt_nums):
		if gt_nums == '1|1':
			return 1.0 
		elif gt_nums == None:
			return 0.0 
		return 0.5

	
	

	def getMatchingVCFRecordFromTrueVCF (self, vcf_record):
		"""
			Returns a matching VCF record from the true vcf file in order to determine
			the true VAFs and whether it is a somatic mutation or not
		"""
		deletion = Deletion(vcf_record)
		 

		candidates = self.true_vcf_reader.fetch(returnChromosome(vcf_record.CHROM), max(0, vcf_record.POS - 250), vcf_record.POS + 250)
		print(len(candidates))
		for candidate in candidates: 
			if isDeletion(candidate):
				continue
			candidate_deletion = Deletion(candidate)
			d = returnDistance(deletion.centerpoints, candidate_deletion.centerpoints)
			if d <= 20 and abs(deletion.length - candidate_deletion.length) <= 10:
				return vcf_record

		return None
			


	def getTruth (self, vcf_record):
		"""	Obtains the ground truth for a VCF record. Input is a vcf record from the 
			aligner; the variant represented by this record may, therefore, not exist. 
			We first match the vcf record with a vcf record from the true vcf file. 
			If there is no match, we claim the vcf record had a frequency of 0.0 in the 
			control and cancer sample.	
			
			Return:
				- call; True, when aligner was (approximately) correct and false otherwise. 
				- h_vaf; VAF of the healthy cells
				- c_vaf; VAF of the cancer cells
				- somatic_mutation; whether is a is a somatic mutation or not
		"""
		if self.truth: 
			h_vaf 		= 0.0  # healthy VAF
			# VAFs of the 4 cancer populations 
			som1_vaf 	= 0.0 
			som2_vaf 	= 0.0 
			som3_vaf 	= 0.0 
			som4_vaf 	= 0.0 

			for call in vcf_record.samples:
				if call.sample=='Control':
					h_vaf = self.determineVAF (call.gt_nums)
				if call.sample=='Som1':
					som1_vaf = self.determineVAF (call.gt_nums)	
				if call.sample=='Som2':
					som2_vaf = self.determineVAF (call.gt_nums)	
				if call.sample=='Som3':
					som3_vaf = self.determineVAF (call.gt_nums)
				if call.sample=='Som4':
					som4_vaf = self.determineVAF (call.gt_nums)
		
			cancer_vaf = (1.0 / 3.0) * som1_vaf + (1.0 / 3.0) * som2_vaf + (1.0 / 6.0) * som3_vaf + (1.0 / 6.0) * som4_vaf 

			somatic_mutation = False
			if h_vaf == 0.0 and cancer_vaf > 0.0:
				somatic_mutation = True

			return True, h_vaf, cancer_vaf, somatic_mutation  

		else:
			


			true_vcf_record = self.getMatchingVCFRecordFromTrueVCF(vcf_record)
			if true_vcf_record == None:
				return False, 0.0, 0.0, False

			h_vaf 		= 0.0  # healthy VAF
			# VAFs of the 4 cancer populations 
			som1_vaf 	= 0.0 
			som2_vaf 	= 0.0 
			som3_vaf 	= 0.0 
			som4_vaf 	= 0.0 

			for call in true_vcf_record.samples:
				if call.sample=='Control':
					h_vaf = self.determineVAF (call.gt_nums)
				if call.sample=='Som1':
					som1_vaf = self.determineVAF (call.gt_nums)	
				if call.sample=='Som2':
					som2_vaf = self.determineVAF (call.gt_nums)	
				if call.sample=='Som3':
					som3_vaf = self.determineVAF (call.gt_nums)
				if call.sample=='Som4':
					som4_vaf = self.determineVAF (call.gt_nums)
		
			cancer_vaf = (1.0 / 3.0) * som1_vaf + (1.0 / 3.0) * som2_vaf + (1.0 / 6.0) * som3_vaf + (1.0 / 6.0) * som4_vaf 

			somatic_mutation = False
			if h_vaf == 0.0 and cancer_vaf > 0.0:
				somatic_mutation = True

			return True, h_vaf, cancer_vaf, somatic_mutation  



	def updateDeletion (self, vcf_record, somatic_mutation, h_mle, c_mle, ci, unique, indel, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments):
		exists, h_true_vaf, c_true_vaf, true_sm = self.getTruth(vcf_record)
		test_case = SMTestCase(True, exists, h_true_vaf, c_true_vaf, true_sm, somatic_mutation, h_mle, c_mle, ci, unique, indel, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments)
		self.testcases.append(test_case)
		self.testcases[-1].print()

	def updateInsertion (self, vcf_record, somatic_mutation, h_mle, c_mle, ci, unique, indel, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments):
		exists, h_true_vaf, c_true_vaf, true_sm = self.getTruth(vcf_record)
		test_case = SMTestCase(False, exists, h_true_vaf, c_true_vaf, true_sm, somatic_mutation, h_mle, c_mle, ci, unique, indel, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments)
		self.testcases.append(test_case)
	
	def printTable2x2 (self, table, title):
		row_total, column_total = [0,0], [0,0]
		n = 0

		for r in [0,1]:
			for c in [0,1]:
				row_total[r] = row_total[r] + table[r][c]
				column_total[c] = column_total[c] + table[r][c]
				n = n + table[r][c]

		recall = None
		precision = None
		if column_total[0] != 0:
			recall = table[0][0] / float(column_total[0])
		if row_total[0] != 0:
			precision = table[0][0] / float(row_total[0])

		print('------------------------TABLE', title, '------------------------')
		print('\t\t\t\t\ttruth')
		print('\t\t|\tsomatic\t\t|\tnot somatic\t|\ttotal')
		print('--------------------------------------------------------------------------------------------')
		print('somatic\t\t|\t', table[0][0], '\t\t|\t', table[0][1], '\t\t|\t', row_total[0])
		print('not somatic\t|\t', table[1][0], '\t\t|\t', table[1][1], '\t\t|\t', row_total[1])	
		print('--------------------------------------------------------------------------------------------')
		print('total\t\t|\t', column_total[0], '\t\t|\t', column_total[1], '\t\t|\t', n)

		print('\nRECALL:\t\t', recall)
		print('PRECISION:\t', precision)	
		print('-------------------------------------------------------------------\n\n')	


	def summarize (self):
		"""Processes the gathered information, generates a summary and outputs it to the command line"""
		for range in self.length_range:
			total = 0.0			# total number of indels in this range 
			n_unique = 0.0			# number of indels for which the MLE was unique -> there was data to say something about their VAF
			
			table = [[0,0],[0,0]]		
			for test_case in self.testcases:
				if test_case.isOfLength(range[0], range[1]):
					total = total + 1
					if test_case.unique: 
						n_unique = n_unique + 1
						if test_case.true_sm and test_case.call:
							table[0][0] = table[0][0] + 1
						elif test_case.true_sm and not test_case.call:
							table[1][0] = table[1][0] + 1
						elif not test_case.true_sm and test_case.call:
							table[0][1] = table[0][1] + 1
						else: 
							table[1][1] = table[1][1] + 1

			print('*************************** LENGTH RANGE', range, '*******************************')
			print('total:\t\t', int(total))
			print('# unique:\t', int(n_unique), '\n')
			self.printTable2x2(table, 'Performance')
			print('*******************************************************************************************************')
	

	def printBoolean (self, boolean, f):
		if boolean:
			print('1\t', end = "", file = f)
		else:
			print('0\t', end = "", file = f)

	def printNumber (self, x, f):
		print(str(x) + '\t', end = "", file = f)

	def printList (self, l, f):
		s = '['
		if len(l) == 0:
			s = '[]'	
		else:
			for i in range(len(l) - 1):
				s = s + str(l[i]) + ', '
			s = s + str(l[-1]) + ']'
		print(s + '\t', end = "", file = f)

	def printToFile (self, output_filename):
		"""
			Stores the collected results into a file. 
		"""		
		f = open(output_filename, 'w')

		print('del\tchr\tpos\tlength\texists\tcorrect\tsomatic\th_vaf\tc_vaf\tunique\tcall\th_mle\tc_mle\tci_start\tci_end\tX_h\tY_h\tX_c\tY_c', file = f) # header

		for t in self.testcases: 
			self.printBoolean(t.is_deletion, f)
			print(str(t.indel.chromosome) + '\t', end = "", file = f)
			self.printNumber(t.position, f)
			self.printNumber(t.length, f)
			self.printBoolean(t.exists, f)
			self.printBoolean(t.correct, f)
			self.printBoolean(t.true_sm, f)
			self.printNumber(t.h_true_vaf, f)
			self.printNumber(t.c_true_vaf, f)
			self.printBoolean(t.unique, f)
			self.printBoolean(t.call, f)
			self.printNumber(t.h_mle, f)
			self.printNumber(t.c_mle, f)
			self.printNumber(t.ci[0], f)
			self.printNumber(t.ci[1], f)
			self.printList(t.Xh, f)
			self.printList(t.Yh, f)
			self.printList(t.Xc, f)	
			self.printList(t.Yc, f)
			print('\n', end="", file = f) 
		f.close()


def main():


	parser = OptionParser(usage=usage)
	
	parser.add_option("-a", action="store", dest="alpha", default=0.0, type=float, 
					help="Level of impurity in disease/cancer sample. (Default = 0.0; no healthy cells present in the disease/cancer sample)")
	parser.add_option("-p", action="store", dest="ploidy", default=None, type=int, 
						help="Ploidy: the number of sets of chromosomes in the nucleus of a cancer cell. (Default = no assumptions made; unit interval is used)")
	parser.add_option("-c", action="store", dest="chromosome", default="", type=str, 
					help="In case one wants to analyze specifically one chromosome, e.g., '-c = 22' would entail that only indels on chr. 22 are consdired. (Default: all chromosomes)")
	parser.add_option("-r", action="store", dest="search_range", default=5000, type=int, 
					help="Range to search for potentially relevant reads (Default: 5000 bp)")
	parser.add_option("--uncertainty_off", action="store_true", dest="uncertainty_off", default=False,
						help="Alignment uncertainties are neglected. (Default = false)")
	parser.add_option("-o", action="store", dest="output_file", default="somatic_mutation_caller_output.txt", 
						help="Output file for the results. Only used for testing. (Default = somatic_mutation_caller_output.txt)")
	parser.add_option("--truth", action="store_true", dest="truth", default=False,
						help="Takes only the true VCF file into account. This is used for testing only since it is not similar to the situation in the lab. (Default = False)")

	(options, args) = parser.parse_args()


	vcf_filename = None
	bam_healthy_filename = None
	bam_cancer_filename = None

	if options.truth:
		if len(args)!=4:
			parser.print_help()
			return 1
		else:
			vcf_filename = "/data1/structvar/vaf-experiments/vcf/cancerclones.vcf.gz"
			bam_healthy_filename = os.path.abspath(args[2])
			bam_cancer_filename = os.path.abspath(args[3])
	else:
		if len(args)!=5:
			parser.print_help()
			return 1
		else: 
			vcf_filename = os.path.abspath(args[2])
			bam_healthy_filename = os.path.abspath(args[3])
			bam_cancer_filename = os.path.abspath(args[4])


	# mean and standard deviation of the normal distribution that the internal segment lengths are assumed to follow when
	# not affected by an indel
	mu 	= float(args[0])
	sigma 	= float(args[1])
	
	# Set the support for the variant allele frequency
	cancer_vaf_range = [] # VAF can vary over the unit interval	
	if options.ploidy is not None: 
		cancer_vaf_range = [x/float(options.ploidy) for x in range(0,options.ploidy+1)] 
	else:
		cancer_vaf_range = None

	sm_caller = SomaticMutationCaller(vcf_filename, bam_healthy_filename, bam_cancer_filename, mu, sigma, alpha = options.alpha, cancer_vaf_range = cancer_vaf_range, search_range = options.search_range, output_filename = options.output_file, no_uncertainty = options.uncertainty_off, truth=options.truth) 
	sm_caller.call() 

if __name__ == '__main__':
	sys.exit(main())

