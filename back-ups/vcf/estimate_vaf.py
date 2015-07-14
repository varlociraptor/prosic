#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
import os
import sys
import vcf

from Variant.Indel import *
from BAM.Alignments import *
from BAM.BAMProcessor import *
from Model.VAFLikelihoodMaximizer import *

__author__ = "Louis Dijkstra"

usage = """%prog [options] <mu> <sigma> <vcf-filename> <bam-filename> 

Estimates the variant allele frequency for the deletions present in 
the VCF file <vcf-file> on the basis of the alignment data in the BAM file
<bam-file>. The results are outputed as an annotated VCF file <output-vcf>.

Have a look at the options '--unit' and '-p' for adjusting the range of the VAF.  

The maximum a posteriori (MAP) estimate of the VAF is added to the INFO field
of the VCF record under the key 'MAP_VAF'. 

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
	<output-vcf> 	annotated VCF file. Contains the variants considered
				extended with the VAF estimate.

NOTE: 	the current implementation ignores indels in the VCF file 
	for which more than one alternative (ALT) are provided. 
	In addition, a uniform prior for the VAF is assumed. 	
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
	parser.add_option("--unit", action="store_true", dest="continuous", default=False,
						help="The variant allele frequency (VAF) can take any value in the unit interval. Overwrites the '-p' option. (Default = False)")
	parser.add_option("-e", action="store", dest="epsilon_a", default=0.0, type=float, 
						help="Value of the parameter epsilon_a (the probability of observing a split when the indel is not present; Default = 0.0)")
	parser.add_option("-E", action="store", dest="epsilon_p", default=0.0, type=float, 
						help="Value of the parameter epsilon_p (the probability of not observing a split when the indel is present; Default = 0.0)")
	parser.add_option("-m", action="store", dest="min_length", default=10, type=int,
				  		help="Minimal length of an indel to be considered. (Default = 10 bp)")
	parser.add_option("-p", action="store", dest="ploidy", default=2, type=int, 
						help="Ploidy: the number of sets of chromosomes in the nucleus of a cell. (Default = 2; diploid)")
	parser.add_option("-r", action="store", dest="search_range", default=5000, type=int, 
						help="Range to search for potentially relevant reads (Default = 5000 bp)")
	parser.add_option("-t", action="store", dest="tolerance", default=0.00000001, type=float, 
						help="Precision threshold used for numerically maximizing the loglikelihood function and computing integrals. (Default = 10^-8)")
	parser.add_option("-x", action="store", dest="autosome", default=None, type=int, 
						help="Processes only this autosome. (Default = all are processed)")
	(options, args) = parser.parse_args()
	
	if (len(args)!=5):
		parser.print_help()
		return 1

	# typing uncertainty for split reads is turned off: 
	if options.typing_uncertainty_off: 
		options.epsilon_a, options.epsilon_p = 0.0, 0.0

	# list of autosomes to be processed 
	list_of_autosomes = range(1,23) # default: all autosomes
	if options.autosome is not None:
		list_of_autosomes = [options.autosome]

	# Set the support for the variant allele frequency
	vaf_range = [] # VAF can vary over the unit interval	
	if not options.continuous: 
		vaf_range = [x/float(options.ploidy) for x in range(0,options.ploidy+1)] 
	
	# mean and standard deviation of the normal distribution that the internal segment lengths are assumed to follow when
	# not affected by an indel
	mu 	= float(args[0])
	sigma 	= float(args[1])
	
	vcf_filename 		= os.path.abspath(args[2])
	bam_filename 		= os.path.abspath(args[3])
	annotated_vcf_filename 	= os.path.abspath(args[4])

	vcf_reader 			= vcf.Reader(open(vcf_filename))
	vcf_writer 			= vcf.Writer(open(annotated_vcf_filename, 'w'), vcf_reader)
	bam_processor 			= BAMProcessor(bam_filename, search_range = options.search_range, primary_alignments_only = options.primary_only) # processes the BAM file
	vaf_likelihood_maximizer	= VAFLikelihoodMaximizer(mu, sigma, epsilon0 = options.epsilon_a, epsilon1 = options.epsilon_p, vaf_range = vaf_range, tolerance = options.tolerance) # deals with the likelihood function 

	# Walk through all records in the given vcf file
	for vcf_record in vcf_reader: 
		if not returnAutosome(vcf_record.CHROM) in list_of_autosomes: # not the autosome we are currently interested in.
			continue
		if len(vcf_record.ALT) != 1: # records with several alternatives are ignored
			continue 
		if not isDeletion(vcf_record): # not a deletion.
			continue 
		if returnIndelLength(vcf_record) < options.min_length: # below the minimal length of an indel.
			continue 

		deletion = Deletion(vcf_record)

		# Obtain the evidence (internal segment based and overlapping alignments): 
		paired_end_alignments, overlapping_alignments = bam_processor.processDeletion(Deletion(vcf_record))
	
		if options.align_uncertainty_off:
			for a in paired_end_alignments: 	a.probability = 1.0 
			for a in overlapping_alignments: 	a.probability = 1.0 

		if options.internal_segment_based_only: # overlapping alignments are turned off 
			for a in overlapping_alignments: 	a.probability = 0.0	

		if options.split_read_alignments_only: # internal segment based evidences are discarded
			for a in paired_end_alignments: 	a.probability = 0.0

		MAP_VAF = vaf_likelihood_maximizer.maximize(vcf_record, paired_end_alignments, overlapping_alignments)
		vcf_record.INFO['MAP_VAF'] = MAP_VAF
		vcf_writer.write_record(vcf_record)

	bam_processor.close()

if __name__ == '__main__':
	sys.exit(main())
