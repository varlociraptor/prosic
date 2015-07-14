#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
import os
import vcf
import pysam

__author__ = "Louis Dijkstra"

from Indel import *
from Alignments import *
from BAMProcessor import *
from SMLikelihoodMaximizer import *

usage = """%prog [options] <mu> <sigma> <vcf-filename> <bam-healthy> <bam-cancer> 

Calls the deletions present in the VCF file as either (1) somatic, (2) germline and 
(3) not present, on the basis of the alignments in the BAM files for the control and the
cancer/disease sample. The healthy cells are assumed to be diploid. 

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

NOTE: the current implementation ignores indels in the vcf file 
		for which more than one alternative (ALT) are provided. 	
"""

# TODO export the results in VCF format
# TODO change to the latest calling method (now confidence interval is still used)
# TODO TODO TODO do this file last! 

def printDeletion(deletion):
	"""Print information about the deletion to standard output"""
	print(deletion.chromosome, '\t', deletion.vcf_record.POS, '\t', deletion.length, '\t', end = "")	

def printAlignments (alignments):
	"""Prints the raw data in tuples: <observation>/<probability>"""
	print('[ ', end = "")
	for a in alignments:
		print(str(a.value) + '/' + str(a.probability) + ' ', end = "")
	print(']', end = "")

def printContinuous(deletion, MLE, unique, call, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments):
	"""Prints the results for a deletion when the support of the cancer VAF is chosen to be continuous"""
	printDeletion(deletion)
	if unique:
		print(MLE, '\t' + call + '\t', end = "")
		printAlignments(h_paired_end_alignments)
		print('\t', end = "")
		printAlignments(h_overlapping_alignments)
		print('\t', end = "")
		printAlignments(c_paired_end_alignments)
		print('\t', end = "")
		printAlignments(c_overlapping_alignments)
		print('') # end line
	else:
		print('Insufficient data; likelihood function has no unique global maximum.')
		

def printDiscrete(deletion, enough_certainty, prob_not_present, prob_somatic, prob_germline, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments):
	"""Prints the results for a deletion when the support of the cancer VAF is chosen to be discrete"""
	printDeletion(deletion)
	if enough_certainty:
		print('1\t', end = "")
	else:
		print('0\t', end = "")
	print(prob_not_present, '\t', end = "")
	print(prob_somatic, '\t', end = "")
	print(prob_germline, '\t', end = "")
	printAlignments(h_paired_end_alignments)
	print('\t', end = "")
	printAlignments(h_overlapping_alignments)
	print('\t', end = "")
	printAlignments(c_paired_end_alignments)
	print('\t', end = "")
	printAlignments(c_overlapping_alignments)
	print('') # end line
	

def determineProbabilities (posterior_distr, cancer_vaf_range):
	"""Returns three probabilities:
		(1) probability that the is not present in both the control and disease sample
		(2) probability that the deletion is somatic
		(3) probability that the the deletion is germline
	   on the basis of the posterior distribution. Only applies when the cancer VAF is assumed discrete."""
	prob_not_present = posterior_distr[0][0] # h_vaf and c_vaf are both zero
	prob_somatic = 0.0
	for i in range(1, len(cancer_vaf_range)): # h_vaf = 0 but c_vaf > 0
		prob_somatic += posterior_distr[0][i]
	return prob_not_present, prob_somatic, 1.0 - prob_not_present - prob_somatic
	

def main():

	parser = OptionParser(usage=usage)
	parser.add_option("--primary-only", action="store_true", dest="primary_only", default=False,
						help="Only primary alignments are taken into account. (Default = False)")
	parser.add_option("--uncertainty-off", action="store_true", dest="uncertainty_off", default=False,
						help="Alignment uncertainty is discared. Every alignment is taken 100% seriously. (Default = False)")
	parser.add_option("-a", action="store", dest="alpha", default=0.0, type=float, 
						help="Level of impurity in disease/cancer sample. (Default = 0.0; no healthy cells present in the disease/cancer sample)")
	parser.add_option("-c", action="store", dest="certainty_threshold", default=0.95, type=float, 
						help="A call is made only when the largest posterior probability is above this threshold. Applies when option '-p' is used. (Default = 0.95)")
	parser.add_option("-m", action="store", dest="min_length", default=10, type=int,
				  		help="Minimal length of an indel to be considered. (Default = 10)")
	parser.add_option("-p", action="store", dest="ploidy", default=None, type=int, 
						help="Ploidy: the number of sets of chromosomes in the nucleus of a cancer cell. When not used, the VAF of the cancer cells can vary over the unit interval.")
	parser.add_option("-r", action="store", dest="search_range", default=5000, type=int, 
						help="Range to search for potentially relevant reads (Default = 5000 bp)")
	parser.add_option("-t", action="store", dest="tolerance", default=0.0000001, type=float, 
						help="Precision threshold used for numerically maximizing the loglikelihood function and the confidence interval (Default = 10^-7)")
	(options, args) = parser.parse_args()
	
	if (len(args)!=5):
		parser.print_help()
		return 1

	# mean and standard deviation of the normal distribution that the internal segment lengths are assumed to follow when
	# not affected by an indel
	mu 	= float(args[0])
	sigma 	= float(args[1])
	
	vcf_filename 		= os.path.abspath(args[2])
	bam_healthy_filename 	= os.path.abspath(args[3])
	bam_cancer_filename 	= os.path.abspath(args[4])
	
	# Set the support for the variant allele frequency
	cancer_vaf_range = None # VAF can vary over the unit interval	
	if options.ploidy is not None: 
		cancer_vaf_range = [x/float(options.ploidy) for x in range(0,options.ploidy+1)] 

	vcf_reader 		= vcf.Reader(open(vcf_filename))
	bam_healthy_processor 	= BAMProcessor(bam_healthy_filename, search_range = options.search_range, primary_alignments_only = options.primary_only) # processes the BAM file of the control/healthy sample
 	bam_cancer_processor 	= BAMProcessor(bam_cancer_filename, search_range = options.search_range, primary_alignments_only = options.primary_only) # processes the BAM file of the cancer sample
 	sm_likelihood_maximizer = SMLikelihoodMaximizer(mu, sigma, alpha = options.alpha, cancer_vaf_range = cancer_vaf_range, tolerance=options.tolerance) # deals with the likelihood function

	# Walk through all records in the VCF file
	for vcf_record in vcf_reader:  
		if len(vcf_record.ALT) != 1: # records with several alternatives are ignored
			continue 
		if not isDeletion(vcf_record): 
			continue 
		if returnIndelLength(vcf_record) < options.min_length:
			continue 
		deletion = Deletion(vcf_record)
		# Obtain the evidence (internal segment based and overlapping alignments): 
		h_paired_end_alignments, h_overlapping_alignments = bam_healthy_processor.processDeletion(deletion)
		c_paired_end_alignments, c_overlapping_alignments = bam_cancer_processor.processDeletion(deletion)

		if options.uncertainty_off:
			for a in h_paired_end_alignments: 	a.probability = 1.0	
			for a in h_overlapping_alignments: 	a.probability = 1.0	
			for a in c_paired_end_alignments: 	a.probability = 1.0	
			for a in c_overlapping_alignments: 	a.probability = 1.0	

		# Call the variant as being somatic or not
		somatic_mutation = False
		if options.ploidy is None: # cancer VAF is continuous on the unit interval
			[mle_h_vaf, mle_c_vaf], max_logl, unique = sm_likelihood_maximizer.computeMLE_continuous (vcf_record, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments)	
			call = "no call"
			if unique:
				if mle_h_vaf == 0.0:
					if sm_likelihood_maximizer.loglikelihood(mle_h_vaf, 0.0) - max_logl + 1.92 < 0:	
						call = "somatic"
					else:
						call = "not present "
				else:
					call = "germline"
			printContinuous(deletion, [mle_h_vaf, mle_c_vaf], unique, call, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments)
		else: 
			posterior_distr = sm_likelihood_maximizer.posteriorDistribution(vcf_record, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments)
			prob_not_present, prob_somatic, prob_germline = determineProbabilities(posterior_distr, cancer_vaf_range)
			if max([prob_not_present, prob_somatic, prob_germline]) >= options.certainty_threshold: # enough certainty to make a claim
				printDiscrete(deletion, True, prob_not_present, prob_somatic, prob_germline, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments)
			else:
				printDiscrete(deletion, False, prob_not_present, prob_somatic, prob_germline, h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments)
						
				
	bam_healthy_processor.close()
	bam_cancer_processor.close()

if __name__ == '__main__':
	sys.exit(main())
