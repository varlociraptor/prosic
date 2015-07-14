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
	p_somatic_mut	- posterior probability of the hyptohesis that this is a somatic mutation 
	p_germline	- posterior probability of the hyptohesis that this is a germline variant
	p_not_present	- posterior probability of the hyptohesis that this variant is not present at all
	CI_lower	- the lower bound of the 95% credible interval
	CI_upper	- the upper bound of the 95% credible interval
	mle_h_vaf	- the maximum likelihood estimate of the healthy VAF (0/0.5/1)
	mle_c_vaf	- the maximum likelihood estimate of the cancer VAF (unit interval)
	c_vaf_mle0 	- MLE of the cancer VAF when the healthy VAF is equal to 0 
	maxlogl0	- maximum loglikelihood when the healthy VAF is restricted to 0 
	c_vaf_mle1 	- MLE of the cancer VAF when the healthy VAF is equal to 1/2 
	maxlogl1	- maximum loglikelihood when the healthy VAF is restricted to 1/2 
	c_vaf_mle2 	- MLE of the cancer VAF when the healthy VAF is equal to 1 
	maxlogl2	- maximum loglikelihood when the healthy VAF is restricted to 1 

	In conclusion, there are always 13 columns that summarize the results for the continuous case.

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

LENGTH_RANGES 		= [[None, None], [10,29], [30,49], [50,69], [70,99], [100, 249], [250, None]]

VCF_TRUTH_FILENAME = "/data1/structvar/vaf-experiments/vcf/cancerclones.vcf" # VCF file that contains the ground truth
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

def returnIndex (chromosome):
	if len(chromosome) > 3 and chromosome[:3] == 'chr':
		chromosome = chromosome[3:]
	if chromosome == 'X' or chromosome == 'x':
		return 22
	if chromosome == 'Y' or chromosome == 'y':
		return 23	
	return int(chromosome)

def main():

	parser = OptionParser(usage=usage)
	parser.add_option("--align-uncertainty-off", action="store_true", dest="align_uncertainty_off", default=False,
						help="Alignment uncertainty is discared. Every alignment is taken 100% seriously.")
	parser.add_option("--internal-segment-based-only", action="store_true", dest="internal_segment_based_only", default=False,
						help="Only internal segment based evidences are taken into account.")
	parser.add_option("--low", action="store_true", dest="low_precision", default=False,
						help="Uses low precision when matching detected indels with the indels in the true VCF file, i.e., the centerpoints must be less or exactly 100 bp away and the length must not differ more than 100 bp. Tihs option overwrites the options '-d' and '-l' (!).")
	parser.add_option("--primary-only", action="store_true", dest="primary_only", default=False,
						help="Only primary alignments are taken into account. (Default = False)")
	parser.add_option("--split-read-alignments-only", action="store_true", dest="split_read_alignments_only", default=False,
						help="Only split read evidences are taken into account.")
	parser.add_option("--truth", action="store_true", dest="truth", default=False,
						help="The indels in the VCF file containing the true genetic variants are called (/data1/structvar/vaf-experiments/vcf/cancerclones.vcf.gz). No VCF file should be provided. This file is normally only used for validating the results on another VCF file.")
	parser.add_option("--typing-uncertainty-off", action="store_true", dest="typing_uncertainty_off", default=False,
						help="Typing uncertainty for the overlapping alignments is turned off. Equivalent to setting epsilon0 and epsilon1 both to 0.0, i.e., '-e 0.0 -E 0.0'.")
	parser.add_option("-a", action="store", dest="alpha", default=0.0, type=float, 
						help="Level of impurity in disease/cancer sample. (Default = 0.0; no healthy cells present in the disease/cancer sample)")
	parser.add_option("-c", action="store", dest="cap", default=0.05, type=float, 
						help="Puts a cap on the alignment probability, e.g., when '-c .05', then the max. alignment probability is 95%. This is used because we believe that no alignment should be believed to much (Default = 0.05)")
	parser.add_option("-d", action="store", dest="distance_threshold", default=10, type=int,
				  		help="Distance between centerpoints of deletions threshold. The distance must be smaller than this in order for the deletions to be deemed similar (in addition, the lengths must be similar, see option '-l'). (Default = 10 bp)")
	#parser.add_option("-ci", action="store", dest="ci_type", default=0.95, type=float, 
	#					help="Type of confidence interval used for making decisions in case the cancer VAF is continuous. E.g., when '-ci 0.95' - an likelihood-ratio based 95% confidence interval is used (Default = 0.95)")
	parser.add_option("-e", action="store", dest="epsilon0", default=0.0, type=float, 
						help="Value of the parameter epsilon0 (Default = 0.0)")
	parser.add_option("-E", action="store", dest="epsilon1", default=0.0, type=float, 
						help="Value of the parameter epsilon1 (Default = 0.0)")
	parser.add_option("-l", action="store", dest="length_threshold", default=20, type=int,
				  		help="Lengths of deletions threshold. The length difference must be smaller than this in order for the deletions to be deemed similar (in addition, the placement of the centerpoints must be similar, see option '-d'). (Default = 20 bp)")
	parser.add_option("-m", action="store", dest="min_length", default=10, type=int,
				  		help="Minimal length of an indel to be considered. (Default = 10 bp)")
	parser.add_option("-p", action="store", dest="ploidy", default=None, type=int, 
						help="Ploidy: the number of sets of chromosomes in the nucleus of a cancer cell. When not used, the VAF of the cancer cells is assumed continuous and can vary over the unit interval.")
	parser.add_option("-r", action="store", dest="search_range", default=5000, type=int, 
						help="Range to search for potentially relevant reads (Default = 5000 bp)")
	parser.add_option("-t", action="store", dest="tolerance", default=0.00000001, type=float, 
						help="Precision threshold used for numerically maximizing the loglikelihood function and the confidence intervals (Default = 10^-8)")
	parser.add_option("-x", action="store", dest="autosome", default=None, type=int, 
						help="Processes only this autosome. (Default = all are processed)")
	
	(options, args) = parser.parse_args()
	
	vcf_filename, bam_healthy_filename, bam_cancer_filename = None, None, None # allocate memory

	if options.truth:
		if (len(args)!=4):
			parser.print_help()
			return 1
		vcf_filename 		= VCF_TRUTH_FILENAME
		bam_healthy_filename 	= os.path.abspath(args[2])
		bam_cancer_filename 	= os.path.abspath(args[3])
	else:
		if (len(args)!=5):
			parser.print_help()
			return 1
		vcf_filename 		= os.path.abspath(args[2])
		bam_healthy_filename 	= os.path.abspath(args[3])
		bam_cancer_filename 	= os.path.abspath(args[4])

		# determine the thresholds used for determining whether an indel is similar to another
		if options.low_precision:
			options.length_threshold, options.distance_threshold = 100, 100

		obtainListOfTrueDeletions(options.min_length)

	if options.typing_uncertainty_off:
		options.epsilon0, options.epsilon1 = 0.0, 0.0 

	coverage = determineCoverage(bam_cancer_filename) # can either be 40 or 80. Is needed for determining the ground truth since it differs from both files

	list_of_autosomes = range(1,23)
	if options.autosome is not None:
		list_of_autosomes = [options.autosome]

	# mean and standard deviation of the normal distribution that the internal segment lengths are assumed to follow when
	# not affected by an indel
	mu 	= float(args[0])
	sigma 	= float(args[1])
	
	# Set the support for the cancer variant allele frequency
	cancer_vaf_range = None # VAF can vary over the unit interval	
	if options.ploidy is not None: 
		cancer_vaf_range = [x/float(options.ploidy) for x in range(0,options.ploidy+1)] 

	vcf_reader 			= vcf.Reader(open(vcf_filename))
	bam_healthy_processor 		= BAMProcessor(bam_healthy_filename, search_range = options.search_range, primary_alignments_only = options.primary_only) # processes the BAM file of the control/healthy sample
 	bam_cancer_processor 		= BAMProcessor(bam_cancer_filename, search_range = options.search_range, primary_alignments_only = options.primary_only) # processes the BAM file of the cancer sample
	vaf_likelihood_maximizer 	= VAFLikelihoodMaximizer(mu, sigma, epsilon0 = options.epsilon0, epsilon1 = options.epsilon1)
 	sm_likelihood_maximizer 	= SMLikelihoodMaximizer(mu, sigma, alpha = options.alpha, epsilon0 = options.epsilon0, epsilon1 = options.epsilon1, cancer_vaf_range = cancer_vaf_range, tolerance=options.tolerance) # deals with the likelihood function

	# Walk through all records in the VCF file
	for vcf_record in vcf_reader:  
		if not returnIndex(vcf_record.CHROM) in list_of_autosomes:
			continue
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
		else: # apply the cap 
			for a in h_paired_end_alignments: 
				if a.probability > 1.0 - options.cap: 
					a.probability = 1.0 - options.cap
			for a in h_overlapping_alignments: 
				if a.probability > 1.0 - options.cap: 
					a.probability = 1.0 - options.cap
			for a in c_paired_end_alignments: 
				if a.probability > 1.0 - options.cap: 
					a.probability = 1.0 - options.cap
			for a in c_overlapping_alignments: 
				if a.probability > 1.0 - options.cap: 
					a.probability = 1.0 - options.cap

		if options.internal_segment_based_only: # overlapping alignments are turned off 
			for a in h_overlapping_alignments: 	a.probability = 0.0	
			for a in c_overlapping_alignments: 	a.probability = 0.0

		if options.split_read_alignments_only: # internal segment based evidences are discarded
			for a in h_paired_end_alignments: 	a.probability = 0.0
			for a in c_paired_end_alignments: 	a.probability = 0.0

		

		# Determine the posterior distribution on the basis of the control BAM file
		posterior_distribution = vaf_likelihood_maximizer.posteriorDistribution(vcf_record, h_paired_end_alignments, h_overlapping_alignments)
		printList(posterior_distribution)

		p_value = None

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

				p_somatic_mutation, p_germline, p_not_present = sm_likelihood_maximizer.determinePosteriorProbabilitiesHypotheses()
				
				mle_h_vaf = 0.0 
				mle_c_vaf = c_vaf_mle0
				if maxlogl1 > maxlogl0 and maxlogl1 > maxlogl2:
					mle_h_vaf = 0.5 
					mle_c_vaf = c_vaf_mle1
				if maxlogl2 > maxlogl0 and maxlogl2 > maxlogl1:
					mle_h_vaf = 1.0
					mle_c_vaf = c_vaf_mle2

				CI = sm_likelihood_maximizer.determineCredibleInterval(mle_h_vaf, mle_c_vaf, alpha = 0.05, n_steps = 100)
				printList([p_somatic_mutation, p_germline, p_not_present, CI[0], CI[1], mle_h_vaf, mle_c_vaf, c_vaf_mle0, maxlogl0, c_vaf_mle1, maxlogl1, c_vaf_mle2, maxlogl2])
			else: # there is no data to go on... Therefore, everything is completely uncertain: 
				print('0\t', end = "") 
				printList([1.0/3.0, 1.0/3.0, 1.0/3.0, 0.0,1.0, None, None, None, '-inf', None, '-inf', None, '-inf'])
				
		else: # cancer VAF is discrete
			posterior_distribution = sm_likelihood_maximizer.posteriorDistribution(vcf_record, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments)
			for item in posterior_distribution: # print the posterior distribution to the screen
				printList(item)

		printRawData(h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments) # Prints the raw data and alignment probabilities to the screen	
		if options.ploidy is None:
			print('\t', p_value, end = "")
		print('') # newline
		
	bam_healthy_processor.close()
	bam_cancer_processor.close()

if __name__ == '__main__':
	sys.exit(main())
