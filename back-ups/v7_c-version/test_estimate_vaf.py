#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
import os
import sys
import vcf

from Indel import *
from Alignments import *
from BAMProcessor import *
from VAFLikelihoodMaximizer import *

__author__ = "Louis Dijkstra"

usage = """%prog [options] <mu> <sigma> <vcf-filename> <bam-filename> 

Used for testing the VAF estimator. The VCF file must contain the true VAF value.

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
	<bam-file>	BAM file containing the 
				alignments. 
				NOTE:File needs to be sorted 
				and indexed. 

NOTE: the current implementation ignores indels in the VCF file for which more than one alternative (ALT) are provided. 	
"""

class TestCase: 
	"""
		Class to store all results for one test case (deletion/insertion)
	"""
	def __init__(self, is_deletion, mle_disc, unique_disc, mle_cont, unique_cont, ci, indel, paired_end_alignments, overlapping_alignments, posterior_distr):
		self.is_deletion 		= is_deletion		
		self.mle_disc 			= mle_disc
		self.unique_disc 		= unique_disc
		self.mle_cont 			= mle_cont
		self.unique_cont 		= unique_cont
		self.ci 			= ci 
		self.true_vaf 			= indel.getTrueVAF() 		
		self.position 			= indel.vcf_record.POS 
		self.length 			= indel.length
		self.n_paired_end_alignments 	= len(paired_end_alignments)
		self.n_overlapping_alignments 	= len(overlapping_alignments)
		self.posterior_distr 		= posterior_distr 

		self.X = []
		self.Y = []

		for a in paired_end_alignments: 
			self.X.append(a.value)

		for a in overlapping_alignments:
			self.Y.append(a.value)

		if unique_disc:
			if mle_disc == self.true_vaf:
				self.correct_disc = True
			else:
				self.correct_disc = False

		# continuous estimate is considered correct when the 95% CI contains the true VAF 
		if ci[0] <= self.true_vaf and ci[1] >= self.true_vaf:
			self.correct_cont = True
		else:
			self.correct_cont = False

	def print(self, min_length = 0):
		"""
			Prints the relevant data of this test case to the command line
		"""
		if self.is_deletion and self.length >= min_length: # In case of a deletion
			print('--------------------------------------------------------------------------------------------------------')
			print('DEL\t|\t', self.true_vaf, '\t|\t', self.position, '\t', self.length) # print info on deletion 
			print('\nAvailable data: ')
			print('\t\t Internal segment lengths :', self.n_paired_end_alignments, '\t', self.X) 
			print('\t\t Overlapping alignments   :', self.n_overlapping_alignments, '\t', self.Y)
			print('\nDiscrete estimate: ')
			if self.unique_disc:
				if self.correct_disc:
					print('\t\tunique\tcorrect\t\t|\t', self.mle_disc, '\t|\t', max(self.posterior_distr), '\t', self.posterior_distr)
				else:
					print('\t\tunique\tincorrect\t|\t', self.mle_disc, '\t|\t', max(self.posterior_distr), '\t', self.posterior_distr)
			else:
				print('\t\tnot unique') 
			print('\nContinuous estimate: ')
			if self.unique_cont:
				if self.correct_cont:
					print('\t\tunique\tcorrect\t\t|\t', self.mle_cont, '\t|\t', self.ci)
				else:
					print('\t\tunique\tincorrect\t|\t', self.mle_cont, '\t|\t', self.ci)
			else:
				print('\t\tnot unique') 
			print('--------------------------------------------------------------------------------------------------------')
			print('\n\n')
					
		elif self.length >= min_length:
			if self.unique:
				if self.correct:
					print('INS\tunique\tV\t', self.position, '\t', self.mle, '\t', self.true_vaf,'\t', self.length, '\t', max(self.posterior_distr), '\t',self.posterior_distr, '\t', self.n_paired_end_alignments, '\t', self.n_overlapping_alignments, '\t', self.X, '\t', self.Y) 
				else:
					print('INS\tunique\tX\t', self.position, '\t', self.mle, '\t', self.true_vaf,'\t', self.length, '\t', max(self.posterior_distr), '\t', self.posterior_distr, '\t', self.n_paired_end_alignments, '\t', self.n_overlapping_alignments, '\t', self.X, '\t', self.Y) 		

			else:
				print('INS\tnot unique\t\t', self.position, '\t', self.mle, '\t', self.true_vaf,'\t', self.length, '\t',max(self.posterior_distr), '\t', self.posterior_distr, '\t', self.n_paired_end_alignments, '\t', self.n_overlapping_alignments) 
					


	def confidentEnough (self, threshold):
		"""
			Returns True when there is a posterior probability larger or equal to the given threshold; otherwise False. 
		"""
		for p in self.posterior_distr:
			if p >= threshold:
				return True
		return False

	def enoughEvidence (self, n):
		"""
			Returns true when the number of relevant reads is equal or exceeds a given n
		"""
		if self.n_paired_end_alignments + self.n_overlapping_alignments >= n:
			return True
		return False

	def isOfLength(self, a, b): 
		"""Determines whether this test case falls within a certain length class"""
		if a <= self.length <= b:
			return True
		return False


		

class TestDataCollector: 
	"""
		Only used for testing. TODO: Remove later
	"""
	def __init__(self):
		self.testcases = [] 
		self.length_range = [[10, 29],[30,49],[50,99],[100,249]]
		self.min_length = 10 # TODO change?
		
	def updateDeletion(self, MLE_disc, unique_disc, MLE_cont, unique_cont, CI, deletion, paired_end_alignments, overlapping_alignments, posterior_distr):	
		self.testcases.append(TestCase(True, MLE_disc, unique_disc, MLE_cont, unique_cont, CI, deletion, paired_end_alignments, overlapping_alignments, posterior_distr))
		self.testcases[-1].print(min_length = self.min_length)

	def updateInsertion(self, MLE_disc, unique_disc, MLE_cont, unique_cont, CI, insertion, paired_end_alignments, overlapping_alignments, posterior_distr):	
		self.testcases.append(TestCase(False, MLE_disc, unique_disc, MLE_cont, unique_cont, CI, insertion, paired_end_alignments, overlapping_alignments, posterior_distr))
		#self.testcases[-1].print(min_length = self.min_length)


	def printTable3x3 (self, table, title):
		row_total = [0,0,0]
		column_total = [0,0,0]
		n = 0 
		# compute marginals
		for r in [0,1,2]:
			for c in [0,1,2]:
				row_total[r] = row_total[r] + table[r][c]
				column_total[c] = column_total[c] + table[r][c]
				n = n + table[r][c]

		correct = 0 
		for i in [0,1,2]:
			correct = correct + table[i][i]

		incorrect = n - correct

		print('------------------------TABLE', title, '------------------------')
		print('\t\t\t\t\t\ttruth')
		print('\t\t\t|\t0.0\t|\t0.5\t|\t1.0\t|\ttotal')
		print('\t--------------------------------------------------------------------------------------------')
		print('\t\t0.0\t|\t',table[0][0], '\t|\t', table[0][1],'\t|\t', table[0][2],'\t|\t', row_total[0])
		print('assign.\t\t0.5\t|\t',table[1][0], '\t|\t', table[1][1],'\t|\t', table[1][2],'\t|\t', row_total[1])
		print('\t\t1.0\t|\t',table[2][0], '\t|\t', table[2][1],'\t|\t', table[2][2],'\t|\t', row_total[2])
		print('\t--------------------------------------------------------------------------------------------')
		print('\t\ttotal\t|\t', column_total[0], '\t|\t', column_total[1], '\t|\t', column_total[2], '\t|\t',  n)
		print('\nCORRECT:\t', correct, ' /', correct / float(n) * 100.0, '%')
		print('INCORRECT:\t', incorrect, ' /', incorrect / float(n) * 100.0, '%\n')

	def summarize(self, confidence_level = 0.0):
		"""Prints a summary of the test results to the command line"""
		for range in self.length_range:
			cont_correct_case = 0.0 
			total = 0.0 
			
			table = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
			for test_case in self.testcases:
				if test_case.isOfLength(range[0], range[1]):
					total = total + 1
					if test_case.unique_disc and test_case.isOfLength(range[0], range[1]) and test_case.confidentEnough(confidence_level):
						i = int(test_case.mle_disc * 2.0)
						j = int(test_case.true_vaf * 2.0)
						table[i][j] = table[i][j] + 1
					if test_case.correct_cont:
						cont_correct_case = cont_correct_case + 1
			title = 'results for indels of lengths between ' + str(range[0]) + ' and ' + str(range[1])
			self.printTable3x3(table, title)
			print('\nContinuous case: in', int(cont_correct_case), 'of', int(total), 'cases, the 95% confidence interval contained the true VAF, i.e.,', cont_correct_case / total * 100, '%\n') 



def printIndel(indel, is_deletion):
	"""Print information about the indel to standard output"""
	if is_deletion: 
		print('DEL\t', end = '')
	else:
		print('INS\t', end = '')
	print(indel.chromosome, '\t', indel.vcf_record.POS, '\t', indel.length, '\t', indel.getTrueVAF(), '\t', end = '') 


def printAlignments (alignments):
	print('[ ', end = "")
	for a in alignments:
		print(str(a.value) + '/' + "{:.2f}".format(a.probability) + ' ', end = "")
	print(']', end = "")

def printContinuous(indel, is_deletion, MLE, unique, CI, paired_end_alignments, overlapping_alignments):
	"""Prints the results for a deletion when the support of the VAF is chosen to be continuous"""
	printIndel(indel, is_deletion)
	if unique:
		print(MLE, '\t', CI[0], '\t', CI[1], '\t', end = "")
		printAlignments(paired_end_alignments)
		print('\t', end = "")
		printAlignments(overlapping_alignments)
		print('') # end line
	else:
		print('Insufficient data; likelihood function has no unique global maximum.')
		

def printDiscrete(indel, is_deletion, MLE, unique, posterior_distr, paired_end_alignments, overlapping_alignments):
	"""Prints the results for a deletion when the support of the VAF is chosen to be discrete"""
	printIndel(indel, is_deletion)
	if unique:
		print(MLE, '\t', end = "")
		print(posterior_distr, '\t', end = "")
		printAlignments(paired_end_alignments)
		print('\t', end = "")
		printAlignments(overlapping_alignments)
		print('') # end line
	else:
		print('Insufficient data; likelihood function has no unique global maximum.')
	
				
def main():

	parser = OptionParser(usage=usage)
	parser.add_option("--unit", action="store_true", dest="continuous", default=False,
						help="The variant allele frequency (VAF) can take any value in the unit interval. Overwrites the '-p' option. (Default = False)")
	parser.add_option("--deletions-only", action="store_true", dest="deletions_only", default=False,
						help="Only deletions are processed.")
	parser.add_option("--insertions-only", action="store_true", dest="insertions_only", default=False,
						help="Only insertions are processed.")
	parser.add_option("--primary-only", action="store_true", dest="primary_only", default=False,
						help="Only primary alignments are taken into account.")
	parser.add_option("--uncertainty-off", action="store_true", dest="uncertainty_off", default=False,
						help="Alignment uncertainty is discared. Every alignment is taken 100% seriously.")
	parser.add_option("-m", action="store", dest="min_length", default=10, type=int,
				  		help="Minimal length of an indel to be considered. (Default = 10)")
	parser.add_option("-p", action="store", dest="ploidy", default=2, type=int, 
						help="Ploidy: the number of sets of chromosomes in the nucleus of a cell. (Default = 2; diploid)")
	parser.add_option("-r", action="store", dest="search_range", default=5000, type=int, 
						help="Range to search for potentially relevant reads (Default = 5000 bp)")
	parser.add_option("-t", action="store", dest="tolerance", default=0.0000001, type=float, 
						help="Precision threshold used for numerically maximizing the loglikelihood function and the confidence interval (Default = 10^-7)")
	(options, args) = parser.parse_args()
	
	if (len(args)!=4):
		parser.print_help()
		return 1

# mean and standard deviation of the normal distribution that the internal segment lengths are assumed to follow when
	# not affected by an indel
	mu 	= float(args[0])
	sigma 	= float(args[1])
	
	vcf_filename = os.path.abspath(args[2])
	bam_filename = os.path.abspath(args[3])
	
	# Set the support for the variant allele frequency
	vaf_range = [] # VAF can vary over the unit interval	
	if not options.continuous: 
		vaf_range = [x/float(options.ploidy) for x in range(0,options.ploidy+1)] 

	vcf_reader = vcf.Reader(open(vcf_filename))
	bam_processor = BAMProcessor(bam_filename, search_range = options.search_range, primary_alignments_only = options.primary_only) # processes the BAM file
	vaf_likelihood_maximizer = VAFLikelihoodMaximizer(mu, sigma, vaf_range = vaf_range, tolerance = options.tolerance) # deals with the likelihood function

	# Walk through all records in the given vcf file
	for vcf_record in vcf_reader: 
		if len(vcf_record.ALT) != 1: # records with several alternatives are ignored
			continue 

		is_deletion 	= isDeletion(vcf_record) 
		is_insertion 	= isInsertion(vcf_record)
		
		if not is_deletion and not is_insertion:
			continue 
		if options.deletions_only and not is_deletion: 
			continue 
		if options.insertions_only and not is_insertion: 
			continue 
		if returnIndelLength(vcf_record) < options.min_length:
			continue 

		paired_end_alignments, overlapping_alignments = None, None # allocate memory
		indel = Indel(vcf_record)

		# Obtain the evidence (internal segment based and overlapping alignments): 
		if is_deletion:
			paired_end_alignments, overlapping_alignments = bam_processor.processDeletion(Deletion(vcf_record))
		else:
			paired_end_alignments, overlapping_alignments = bam_processor.processInsertion(Insertion(vcf_record))

		if options.uncertainty_off:
			for a in paired_end_alignments: a.probability = 1.0 
			for a in overlapping_alignments: a.probability = 1.0 

		# Estimate the VAF and print the result to the standard output
		if options.continuous:
			MLE, unique, CI = vaf_likelihood_maximizer.computeMLE_95CI (vcf_record, paired_end_alignments, overlapping_alignments)
			printContinuous (indel, is_deletion, MLE, unique, CI, paired_end_alignments, overlapping_alignments)
		else:
			MLE, unique = vaf_likelihood_maximizer.computeMLE_discrete (vcf_record, paired_end_alignments, overlapping_alignments)
			posterior_distr = vaf_likelihood_maximizer.posteriorDistribution(vcf_record, paired_end_alignments, overlapping_alignments)
			printDiscrete (indel, is_deletion, MLE, unique, posterior_distr, paired_end_alignments, overlapping_alignments)

	bam_processor.close()

if __name__ == '__main__':
	sys.exit(main())

