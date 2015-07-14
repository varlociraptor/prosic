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

Estimates the variant allele frequency for the indels (deletion/insertions) 
present in the VCF file and outputs the results.

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
	

# TODO export the results in VCF format

def printIndel(indel, is_deletion):
	"""Print information about the indel to standard output"""
	if is_deletion: 
		print('DEL\t', end = '')
	else:
		print('INS\t', end = '')
	print(indel.chromosome, '\t', indel.vcf_record.POS, '\t', indel.length, '\t', end = '') 

def printDeletion(deletion):
	"""Print information about the deletion to standard output"""
	print(deletion.chromosome, '\t', deletion.vcf_record.POS, '\t', deletion.length, '\t', end = "")	

def printInsertion(insertion):
	"""Print information about the insertion to standard output"""
	print(insertion.chromosome, '\t', insertion.vcf_record.POS, '\t', insertion.length, '\t', end = "")	

def printAlignments (alignments):
	"""Prints the raw data in tuples: <observation>/<probability>"""
	print('[ ', end = "")
	for a in alignments:
		print(str(a.value) + '/' + str(a.probability) + ' ', end = "")
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
