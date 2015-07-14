#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
import os
import sys
import vcf
import pysam
import numpy

from Indel import *
from Alignments import *
from BAMProcessor import *
from VAFLikelihoodMaximizer import * 
from SMLikelihoodMaximizer import *

__author__ = "Louis Dijkstra"

usage = """
When the ground-truth VCF file (/data1/structvar/vaf-experiments/vcf/cancerclones.vcf.gz) is used, call: 

%prog [options] --truth <mu> <sigma> <bam-healthy> <bam-cancer>

Otherwise: 

%prog [options] <mu> <sigma> <vcf-file> <bam-healthy> <bam-cancer> 

Used for testing our capabilities to classify an indel as being (1) not present, (2) somatic and (3) germline on the 
basis of two BAM files; one obtained from a control/healthy sample assumed to be diploid and one based on a sample of the 
cancer tumor. 

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

NOTE: the current implementation ignores indels in the VCF file for which more than one alternative (ALT) are provided. 	

-------------
OUTPUT FORMAT
-------------
Every line represents a deletion of at least a certain minimal length (see option '-m'). The data points are seperated by tabs. 

The first five columns provide info on the deletion: 
	
chromsome 		- the chromosome on which the deletion (is thought to) occur(s)
position   		- reference position as given in the VCF file 
length			- the length of the deletion 
true healthy vaf	- the true variant allele frequency of healthy cells (can either be 0.0, 0.5 or 1.0)
true cancer vaf		- the true variant allele frequency of the cancer cells 

The next three columns provide the posterior distribution of the healthy VAF on the basis of the Control BAM file only (!): 

posterior h_vaf = 0	- posterior probability of being homozygous for NOT having the deletion 
posterior h_vaf = 0.5	- posterior probability heterozygosity 
posterior h_vaf = 1.0	- posterior probability of being homozygous for having the deletion 

So far, the output files are always the same. For the rest it depends on whether the cancer VAF is assumed to vary continuously over 
the unit interval (as is common) or whether one wants to impose a (abitrary) discretization (see option '-p'). In the continuous case: 

Cancer VAF is continuous:
	For every healthy VAF, i.e., 0.0, 0.5 and 1.0, we output:
		c_vaf_mle	- the MLE estimate for the cancer VAF
		maxlogl		- the maximum loglikelihood 
		CI_95_0_down	- the lower bound of the 95% confidence interval
		CI_95_0_up	- the upper bound of the 95% confidence interval
		CI_90_0_down	- the lower bound of the 90% confidence interval
		CI_90_0_up	- the upper bound of the 90% confidence interval
	In conclusion, there are always 18 columns that summarize the results for the continuous case.

Cancer VAF is assumed discrete (option '-p' is used):
	All posterior probabilities are printed. Suppose the cancer VAF can vary over the values c_0, c_1, ..., c_p. 
	First column is the posterior probability of (h_vaf, c_vaf) = (0,c_0), second column (0.5, c_0) and third: (1.0, c_0).
	We then continue: (0.0, c_1), (0.5, c_1), (1.0,c_1) and so forth until c_p. 
	
The last four columns denote the raw data (an array in brackets []) 

	- internal segments control sample
	- overlapping alignment (0/1) control sample
	- internal segments cancer sample
	- overlapping alignment (0/1) cancer sample

Every data point is presented as
	<value>/<alignment probability>
"""

VCF_TRUTH_FILENAME = "/data1/structvar/vaf-experiments/vcf/cancerclones.vcf.gz" # VCF file that contains the ground truth
VCF_TRUTH_READER = vcf.Reader(open(VCF_TRUTH_FILENAME)) # used for reading the VCF file 

list_true_deletions = [[] for i in range(24)] # 24 empty lists for the 22 autosomes and the two sex-chromosomes X and Y   

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
		if not isDeletion(true_vcf_record): 
			continue 
		if returnIndelLength(true_vcf_record) < min_length:
			continue 
		list_true_deletions[returnIndex(true_vcf_record.CHROM)].append(true_vcf_record)

def retrieveMatchingTrueVCFRecord (vcf_record, distance_threshold, length_threshold):
	"""Returns a matching vcf record from the true VCF file"""
	index = returnIndex(vcf_record.CHROM)
	deletion = Deletion(vcf_record)
	candidates = [] # list with potentially matching candidates
	for true_vcf_record in list_true_deletions[index]:
		if abs(true_vcf_record.POS - vcf_record.POS) > 250:
			continue
		true_deletion = Deletion(true_vcf_record)
		if returnMinimumDifference(true_deletion.centerpoints, deletion.centerpoints) <= distance_threshold and abs(deletion.length - true_deletion.length) <= length_threshold:
			return true_vcf_record
	return None
		

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

def obtainTruth(vcf_record, true_vcf_file_used, coverage, distance_threshold, length_threshold):
	"""Returns the true VAFs of the healthy and cancer cells. The first argument (boolean) denotes whether the true VCF file was used
	   for determining the indels, since then it is easy. When false, we first need to locate the corresponding truely present genetic variant and 
	   retrieve the truth. The coverage determines with what frequency certain cancer populations are present."""
	if true_vcf_file_used: 
		return returnTrueVAFs(coverage, vcf_record)
	else: # other VCF than the truth is used
		true_vcf_record = retrieveMatchingTrueVCFRecord (vcf_record, distance_threshold, length_threshold) # retrieve true_vcf_file	
		if true_vcf_record is None:
			return 0.0, 0.0 # deletion does not appear 
		return returnTrueVAFs(coverage, true_vcf_record)

def determineCoverage(bam_cancer_filename):
	"""Determines the coverage used for the simulated cancer files by checking the filename."""
	if 'Cancer40' in bam_cancer_filename:
		return 40 
	elif 'Cancer80' in bam_cancer_filename:
		return 80
	else:
		print('ERROR: Unknown BAM files. Ground truth is unkown and, therefore, no validation is possible. Check given BAM file.')
		sys.exit()

def printList(l):
	"""Outputs a list to standard output. Every element is seperated by a tab."""
	for item in l:
		print(item, '\t', end = "")	

def printAlignmentList (alignment_list):
	"""Outputs a list of alignment data in the form <value>/<alignment probability>.""" 
	print('[ ', end = "")
	for a in alignment_list:
		print(str(a.value) + '/' + str(a.probability) + ' ', end = "")
	print(']', end = "")

def printRawData(h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments):
	"""Prints the raw alignment data to standard output. Every observation is printed as <value>/<alignment probability>."""
	printAlignmentList(h_paired_end_alignments) 
	print('\t', end = "")
	printAlignmentList(h_overlapping_alignments) 
	print('\t', end = "")
	printAlignmentList(c_paired_end_alignments) 
	print('\t', end = "")
	printAlignmentList(c_overlapping_alignments) 


def decodeData(data_vector):
	"""Decodes the data_vector [ <value_1/align_prob_1> <value_2/align_prob_2> ... <value_n/align_prob_n>]
	   Returns a list with the values and a list with the associated alignment probabilities. """
	data_vector = data_vector[1:-1] # throw away brackets []
	values 		= [] # list with data values  
	align_probs 	= [] # list with associated alignment probabilities
	for data_point in data_vector.split(' '):
		if data_point == "" or data_point == '[' or data_point == ']':
			continue
		d = data_point.split('/')
		values.append(d[0])
		align_probs.append(d[1])
	return values, align_probs

def main():

	parser = OptionParser(usage=usage)
	parser.add_option("--align-uncertainty-off", action="store_true", dest="align_uncertainty_off", default=False,
						help="Alignment uncertainty is discared. Every alignment is taken 100% seriously.")
	parser.add_option("--cap-off", action="store_true", dest="cap_off", default=False,
						help="The original alignment probability of the aligner are used. Normally they are capped at 95%.")
	parser.add_option("--typing-uncertainty-off", action="store_true", dest="typing_uncertainty_off", default=False,
						help="Typing uncertainty for the overlapping alignments is turned off. Equivalent to setting epsilon0 and epsilon1 both to 0.0, i.e., '-e 0.0 -E 0.0'.")
	parser.add_option("-a", action="store", dest="alpha", default=0.0, type=float, 
						help="Level of impurity in disease/cancer sample. (Default = 0.0; no healthy cells present in the disease/cancer sample)")
	parser.add_option("-e", action="store", dest="epsilon0", default=0.0, type=float, 
						help="Value of the parameter epsilon0 (Default = 0.0)")
	parser.add_option("-E", action="store", dest="epsilon1", default=0.0, type=float, 
						help="Value of the parameter epsilon1 (Default = 0.0)")
	parser.add_option("-p", action="store", dest="ploidy", default=None, type=int, 
						help="Ploidy: the number of sets of chromosomes in the nucleus of a cancer cell. When not used, the VAF of the cancer cells is assumed continuous and can vary over the unit interval.")
	parser.add_option("-t", action="store", dest="tolerance", default=0.00000001, type=float, 
						help="Precision threshold used for numerically maximizing the loglikelihood function and the confidence intervals (Default = 10^-8)")
	(options, args) = parser.parse_args()
	
	if (len(args)!=3):
		parser.print_help()
		return 1

	# mean and standard deviation of the normal distribution that the internal segment lengths are assumed to follow when
	# not affected by an indel
	mu 			= float(args[0])
	sigma 			= float(args[1])
	preprocessed_file 	= os.path.abspath(args[2])

	if options.typing_uncertainty_off:
		options.epsilon0, options.epsilon1 = 0.0, 0.0 

	# Set the support for the cancer variant allele frequency
	cancer_vaf_range = None # VAF can vary over the unit interval	
	if options.ploidy is not None: 
		cancer_vaf_range = [x/float(options.ploidy) for x in range(0,options.ploidy+1)] 

	vaf_likelihood_maximizer 	= VAFLikelihoodMaximizer(mu, sigma, epsilon0 = options.epsilon0, epsilon1 = options.epsilon1)
 	sm_likelihood_maximizer 	= SMLikelihoodMaximizer(mu, sigma, alpha = options.alpha, epsilon0 = options.epsilon0, epsilon1 = options.epsilon1, cancer_vaf_range = cancer_vaf_range, tolerance=options.tolerance) # deals with the likelihood function

	# Walk through all indels in the preprocessed file 
	for line in preprocessed_file:
		values = line.split('\t') # all values are stored here
		is_deletion = False
		if values[0] == 'DEL' or values[0] == 'DEL ':
			is_deletion = True
		chromosome		= values[1][3:]
		position 		= int(values[2])
		length 		= int(values[3])
		true_h_vaf	= float(values[4])
		true_c_vaf 	= float(values[5])
		Xh, Xh_p	= decodeData(values[6])
		Yh, Yh_p	= decodeData(values[7])
		Xc, Xc_p	= decodeData(values[8])
		Yc, Yc_p	= decodeData(values[9])

		if options.align_uncertainty_off:
			for a in Xh_p: 	a = 1.0	
			for a in Yh_p: 	a = 1.0	
			for a in Xc_p: 	a = 1.0		
			for a in Yc_p: 	a = 1.0	
		
		delta = length
		if not is_deletion:
			delta *= -1.0


		h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments = [], [], [], []
		for i in range(len(Xh)):
			h_paired_end_alignments.append(DummyAlignment(Xh[i], Xh_p[i]))
		for i in range(len(Yh)):
			h_overlapping_alignments.append(DummyAlignment(Yh[i], Yh_p[i]))
		for i in range(len(Xc)):
			c_paired_end_alignments.append(DummyAlignment(Xc[i], Xc_p[i]))
		for i in range(len(Yc)):
			c_overlapping_alignments.append(DummyAlignment(Yc[i], Yc_p[i]))
			
		vaf_likelihood_maximizer.set(h_paired_end_alignments, h_overlapping_alignments, delta)
		sm_likelihood_maximizer.set(h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments, delta)	

		

	# Walk through all records in the VCF file
	for vcf_record in vcf_reader:  
		if len(vcf_record.ALT) != 1: # records with several alternatives are ignored
			continue 
		if not isDeletion(vcf_record): 
			continue 
		if returnIndelLength(vcf_record) < options.min_length:
			continue 

		deletion = Deletion(vcf_record)
		deletion.print() # Prints the relevant information for the deletion 

		true_h_vaf, true_c_vaf = obtainTruth(vcf_record, options.truth, coverage, options.distance_threshold, options.length_threshold) # obtains the true VAF of the healthy and cancer cells
		print(true_h_vaf, '\t', true_c_vaf, '\t', end = "") # Print the truth 

		# Obtain the evidence (internal segment based and overlapping alignments): 
		h_paired_end_alignments, h_overlapping_alignments = bam_healthy_processor.processDeletion(deletion)
		c_paired_end_alignments, c_overlapping_alignments = bam_cancer_processor.processDeletion(deletion)

		if options.align_uncertainty_off:
			for a in h_paired_end_alignments: 	a.probability = 1.0	
			for a in h_overlapping_alignments: 	a.probability = 1.0	
			for a in c_paired_end_alignments: 	a.probability = 1.0	
			for a in c_overlapping_alignments: 	a.probability = 1.0

		# Determine the posterior distribution on the basis of the control BAM file
		posterior_distribution = vaf_likelihood_maximizer.posteriorDistribution(vcf_record, h_paired_end_alignments, h_overlapping_alignments)
		printList(posterior_distribution)

		if options.ploidy is None: # cancer VAF is continuous on the unit interval 
			sm_likelihood_maximizer.set(h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments, getDelta(vcf_record))
			if sm_likelihood_maximizer.globalMaximumExists():
				print('1\t', end = "")
				# various MLE estimates of the cancer VAF (c_vaf) for various values of h_vaf
				c_vaf_mle0 = fminbound(sm_likelihood_maximizer.help_loglikelihood0, 0.0, 1.0, [], xtol=options.tolerance) # MLE for h_vaf = 0.0		
				c_vaf_mle1 = fminbound(sm_likelihood_maximizer.help_loglikelihood1, 0.0, 1.0, [], xtol=options.tolerance) # MLE for h_vaf = 0.5	
				c_vaf_mle2 = fminbound(sm_likelihood_maximizer.help_loglikelihood2, 0.0, 1.0, [], xtol=options.tolerance) # MLE for h_vaf = 1.0	
				# maxima loglikelihoods for various values of h_vaf
				maxlogl0 = sm_likelihood_maximizer.loglikelihood(0.0, c_vaf_mle0) 
				maxlogl1 = sm_likelihood_maximizer.loglikelihood(0.5, c_vaf_mle1)
				maxlogl2 = sm_likelihood_maximizer.loglikelihood(1.0, c_vaf_mle2)
				# determine the 95% confidence intervals for the c_vaf when h_vaf varies
				CI_95_0 = sm_likelihood_maximizer.computeCI(0.0, c_vaf_mle0, maxlogl0, .95)
				CI_95_1 = sm_likelihood_maximizer.computeCI(0.5, c_vaf_mle1, maxlogl1, .95)
				CI_95_2 = sm_likelihood_maximizer.computeCI(1.0, c_vaf_mle2, maxlogl2, .95)
				# determine the 90% confidence intervals for the c_vaf when h_vaf varies
				CI_90_0 = sm_likelihood_maximizer.computeCI(0.0, c_vaf_mle0, maxlogl0, .90)
				CI_90_1 = sm_likelihood_maximizer.computeCI(0.5, c_vaf_mle1, maxlogl1, .90)
				CI_90_2 = sm_likelihood_maximizer.computeCI(1.0, c_vaf_mle2, maxlogl2, .90)
				
				printList([c_vaf_mle0, maxlogl0, CI_95_0[0], CI_95_0[1], CI_90_0[0], CI_90_0[1]])
				printList([c_vaf_mle1, maxlogl1, CI_95_1[0], CI_95_1[1], CI_90_1[0], CI_90_1[1]])
				printList([c_vaf_mle2, maxlogl2, CI_95_2[0], CI_95_2[1], CI_90_2[0], CI_90_2[1]])
			else: # there is no data to go on... Therefore, everything is completely uncertain: 
				print('0\t', end = "")
				printList([None, '-inf', 0, 1, 0, 1])
				printList([None, '-inf', 0, 1, 0, 1])
				printList([None, '-inf', 0, 1, 0, 1])
				
		else: # cancer VAF is discrete
			posterior_distribution = sm_likelihood_maximizer.posteriorDistribution(vcf_record, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments)
			for item in posterior_distribution: # print the posterior distribution to the screen
				printList(item)


		printRawData(h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments) # Prints the raw data and alignment probabilities to the screen	
		print('') # newline
		
	bam_healthy_processor.close()
	bam_cancer_processor.close()

if __name__ == '__main__':
	sys.exit(main())
