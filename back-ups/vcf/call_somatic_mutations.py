#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
import os
import vcf
import pysam

__author__ = "Louis Dijkstra"

from Variant.Indel import *
from BAM.Alignments import *
from BAM.BAMProcessor import *
from Model.SMLikelihoodMaximizer import *

usage = """%prog [options] <mu> <sigma> <vcf-file> <bam-healthy> <bam-cancer>  <output-vcf>

Calls the deletions present in the VCF file <vcf-file> as either (1) somatic, (2) germline 
or (3) not present, on the basis of the alignment data in the BAM files for the 
healthy/control sample <bam-healthy> and the cancer/disease sample <bam-cancer>. In addition,
a maximum a posteriori (MAP) estimate of the variant allele frequencies (VAFs) of the healthy and 
the cancer cells are given. The results are outputed as an annotated VCF file <output-vcf>.
 
The results are added to the INFO field of the VCF records, i.e.,: 

'MAP_CANCER_VAF' 	- the MAP estimate of the cancer cells (varies over the unit interval)
'MAP_HEALTHY_VAF' 	- the MAP estimates of the healthy cells. Value can either be 0.0, 0.5 
			  and 1.0, corresponding to 'absent', 'heterozygous' and 'homozygous'. 
'CALL'			- either 'somatic', 'germline' or 'not present'; the hypthesis with the 
			  highest posterior probability. 
'POSTERIOR_PROB' 	- the posterior probability that belong to the call made. 

In the case that there is no data to make a claim, CALL is equal to 'unknown'. Every other field is 
filled with 'none'. 


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
	<output-vcf> 	annotated VCF file. Contains the variants considered
			from <vcf-file> extended with the made call.

NOTE: 	The current implementation ignores indels in the vcf file 
	for which more than one alternative (ALT) are provided.
	In addition, we assume a uniform prior for the VAFs. Healthy cells
	are assumed to be diploid.  	
"""

def returnAutosome (autosome):
	"""Returns an integer denoting the chromosome."""
	if len(autosome) > 3 and autosome[:3] == 'chr':
		autosome = autosome[3:]
	if autosome == 'X' or autosome == 'x':
		return 23
	if autosome == 'Y' or autosome == 'y':
		return 24	
	return int(autosome)

def determineCall (p_somatic, p_germline, p_not_present):  
	"""Determines on the basis of the posterior probabilities the call, i.e., 'somatic', 'germline' or 'not_present'."""
	if p_somatic > p_germline and p_somatic > p_not_present:
		return 'somatic', p_somatic
	if p_germline > p_somatic and p_germline > p_not_present:
		return 'germline', p_germline
	return 'not_present', p_not_present

def main():

	parser = OptionParser(usage=usage)
	parser.add_option("--align-uncertainty-off", action="store_true", dest="align_uncertainty_off", default=False,
						help="Alignment uncertainty is discared. Every alignment is taken for granted.")
	parser.add_option("--internal-segment-based-only", action="store_true", dest="internal_segment_based_only", default=False,
						help="Only internal segment based evidences are taken into account.")
	parser.add_option("--primary-only", action="store_true", dest="primary_only", default=False,
						help="Only primary alignments are taken into account. (Default = False; also secondary etc. alignments are considered.)")
	parser.add_option("--split-read-alignments-only", action="store_true", dest="split_read_alignments_only", default=False,
						help="Only split-read evidences are taken into account.")
	parser.add_option("--typing-uncertainty-off", action="store_true", dest="typing_uncertainty_off", default=False,
						help="Typing uncertainty for the overlapping alignments is turned off. Equivalent to setting epsilon_a and epsilon_p both to 0.0, i.e., '-e 0.0 -E 0.0'.")
	parser.add_option("-a", action="store", dest="alpha", default=0.0, type=float, 
						help="Level of impurity in disease/cancer sample. (Default = 0.0; no healthy cells present in the disease/cancer sample)")
	parser.add_option("-e", action="store", dest="epsilon_a", default=0.0, type=float, 
						help="Value of the parameter epsilon_a (the probability of observing a split when the indel is not present; Default = 0.0)")
	parser.add_option("-E", action="store", dest="epsilon_p", default=0.0, type=float, 
						help="Value of the parameter epsilon_p (the probability of not observing a split when the indel is present; Default = 0.0)")
	parser.add_option("-m", action="store", dest="min_length", default=10, type=int,
				  		help="Minimal length of an indel to be considered. (Default = 10 bp)")
	parser.add_option("-r", action="store", dest="search_range", default=5000, type=int, 
						help="Range to search for potentially relevant reads (Default = 5000 bp)")
	parser.add_option("-t", action="store", dest="tolerance", default=0.00000001, type=float, 
						help="Precision threshold used for numerically maximizing the loglikelihood function and computing integrals. (Default = 10^-8)")
	parser.add_option("-x", action="store", dest="autosome", default=None, type=int, 
						help="Processes only this autosome. (Default = all are processed)")
	(options, args) = parser.parse_args()
	
	if (len(args)!=6):
		parser.print_help()
		return 1

	# typing uncertainty for split reads is turned off: 
	if options.typing_uncertainty_off: 
		options.epsilon_a, options.epsilon_p = 0.0, 0.0 

	# list of autosomes to be processed 
	list_of_autosomes = range(1,23) # default: all autosomes
	if options.autosome is not None:
		list_of_autosomes = [options.autosome]

	# mean and standard deviation of the normal distribution that the internal segment lengths are assumed to follow when
	# not affected by an indel
	mu 	= float(args[0])
	sigma 	= float(args[1])
	
	vcf_filename 		= os.path.abspath(args[2])
	bam_healthy_filename 	= os.path.abspath(args[3])
	bam_cancer_filename 	= os.path.abspath(args[4])
	annotated_vcf_filename 	= os.path.abspath(args[5])

	vcf_reader 			= vcf.Reader(open(vcf_filename))
	vcf_writer 			= vcf.Writer(open(annotated_vcf_filename, 'w'), vcf_reader)
	bam_healthy_processor 		= BAMProcessor(bam_healthy_filename, search_range = options.search_range, primary_alignments_only = options.primary_only) # processes the BAM file of the control/healthy sample
 	bam_cancer_processor 		= BAMProcessor(bam_cancer_filename, search_range = options.search_range, primary_alignments_only = options.primary_only) # processes the BAM file of the cancer sample
 	sm_likelihood_maximizer 	= SMLikelihoodMaximizer(mu, sigma, alpha = options.alpha, epsilon0 = options.epsilon_a, epsilon1 = options.epsilon_p, cancer_vaf_range = None, tolerance=options.tolerance) # deals with the likelihood function


	# Walk through all records in the VCF file
	for vcf_record in vcf_reader:  
		if not returnAutosome(vcf_record.CHROM) in list_of_autosomes: # not the autosome we are currently interested in.
			continue
		if len(vcf_record.ALT) != 1: # records with several alternatives are ignored.
			continue 
		if not isDeletion(vcf_record): # not a deletion.
			continue 
		if returnIndelLength(vcf_record) < options.min_length: # below the minimal length of an indel.
			continue 

		deletion = Deletion(vcf_record)

		# Obtain the evidence (internal segment based and overlapping alignments): 
		h_paired_end_alignments, h_overlapping_alignments = bam_healthy_processor.processDeletion(deletion)
		c_paired_end_alignments, c_overlapping_alignments = bam_cancer_processor.processDeletion(deletion)

		if options.align_uncertainty_off: # reset the alignment probabilities to 1
			for a in h_paired_end_alignments: 	a.probability = 1.0	
			for a in h_overlapping_alignments: 	a.probability = 1.0	
			for a in c_paired_end_alignments: 	a.probability = 1.0	
			for a in c_overlapping_alignments: 	a.probability = 1.0	

		if options.internal_segment_based_only: # overlapping alignments are turned off 
			for a in h_overlapping_alignments: 	a.probability = 0.0	
			for a in c_overlapping_alignments: 	a.probability = 0.0

		if options.split_read_alignments_only: # internal segment based evidences are discarded
			for a in h_paired_end_alignments: 	a.probability = 0.0
			for a in c_paired_end_alignments: 	a.probability = 0.0
			
		# Call the variant as being somatic/germline or not present
		sm_likelihood_maximizer.set(h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments, getDelta(vcf_record)) # transfers the relevant data to the module dealing with the likelihood function
		
		if sm_likelihood_maximizer.globalMaximumExists(): # the likelihood function is not uniform
			# determine the MAP estimates for the healthy and cancer VAFs
			c_vaf_mle0 = fminbound(sm_likelihood_maximizer.help_loglikelihood0, 0.0, 1.0, [], xtol=options.tolerance) # healthy VAF := 0.0 
			c_vaf_mle1 = fminbound(sm_likelihood_maximizer.help_loglikelihood1, 0.0, 1.0, [], xtol=options.tolerance) # healthy VAF := 0.5
			c_vaf_mle2 = fminbound(sm_likelihood_maximizer.help_loglikelihood2, 0.0, 1.0, [], xtol=options.tolerance) # healthy VAF := 1.0 
			# the loglikelihoods that belong to these maxima:
			maxlogl0 = sm_likelihood_maximizer.loglikelihood(0.0, c_vaf_mle0) 
			maxlogl1 = sm_likelihood_maximizer.loglikelihood(0.5, c_vaf_mle1)
			maxlogl2 = sm_likelihood_maximizer.loglikelihood(1.0, c_vaf_mle2)

			map_h_vaf, map_c_vaf = 0.0, c_vaf_mle0 # initialize the MAP estimates
			if maxlogl1 > maxlogl0 and maxlogl1 > maxlogl2:
				map_h_vaf = 0.5 
				map_c_vaf = c_vaf_mle1
			if maxlogl2 > maxlogl0 and maxlogl2 > maxlogl1:
				map_h_vaf = 1.0
				map_c_vaf = c_vaf_mle2

			# determine the posterior distribution for the hypotheses somatic/germline/not present. 
			p_somatic, p_germline, p_not_present = sm_likelihood_maximizer.determinePosteriorProbabilitiesHypotheses()
			call, posterior_prob = determineCall (p_somatic, p_germline, p_not_present)

			vcf_record.INFO['CALL'] 		= call 
			vcf_record.INFO['POSTERIOR_PROB'] 	= posterior_prob 
			vcf_record.INFO['MAP_HEALTHY_VAF'] 	= map_h_vaf
			vcf_record.INFO['MAP_CANCER_VAF'] 	= map_c_vaf
		else:
			# no statement can be made about the deletion; output an unannotated VCF record. 
			vcf_record.INFO['CALL'] 		= 'unknown'
			vcf_record.INFO['POSTERIOR_PROB'] 	= 'none'
			vcf_record.INFO['MAP_HEALTHY_VAF'] 	= 'none'
			vcf_record.INFO['MAP_CANCER_VAF'] 	= 'none'

		print(vcf_record.INFO) # TODO remove
		vcf_writer.write_record(vcf_record) 	# writes the annotated VCF record to the output VCF file	
	
	vcf_writer.close()
	bam_healthy_processor.close()
	bam_cancer_processor.close()

if __name__ == '__main__':
	sys.exit(main())
