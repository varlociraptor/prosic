#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
import os
import sys
import vcf
import pysam
import numpy as np

__author__ = "Louis Dijkstra"

from Variant.Indel import *
from BAM.Alignments import *
from BAM.BAMProcessor import *

usage = """%prog [options] <vcf-file> <bam-healthy> <bam-cancer> 

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

Outputs the observations for every line in the VCF file <vcf-file> 
obtained from the  alignment data of in the BAM files for the 
healthy/control sample <bam-healthy> and the cancer/disease sample 
<bam-cancer>. For every indel in the VCF file, the program outputs
9 lines. The first contains information on the variant, formatted 
as follows:

	<type>\t<chr>\t<pos>\t<length>

<type> is '+' for an insertion, '-' for a deletion. 
<chr> the autosome on which the variant lies
<pos> position as given in the VCF file
<length> lenght of the indel

The next four lines represent the data of the healthy sample. The 
later four are for the cancer sample. The first line represents in
both cases the observed insert sizes. The second the associated 
alignment probabilities. The third the split read observations (0/1).
The fourth contains probabilities again. 
	
NOTE: 	The current implementation ignores indels in the vcf file 
	for which more than one alternative (ALT) are provided.	
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
	parser.add_option("--deletions-only", action="store_true", dest="deletions_only", default=False, 
				  		help="Only deletions are stored and insertions are discarded.")
	parser.add_option("--insertions-only", action="store_true", dest="insertions_only", default=False, 
				  		help="Only insertions are stored and deletions are discarded.")
	parser.add_option("--primary-only", action="store_true", dest="primary_only", default=False,
						help="Only primary alignments are taken into account. (Default = False; also secondary etc. alignments are considered.)")
	parser.add_option("-m", action="store", dest="min_length", default=10, type=int,
				  		help="Minimal length of an indel to be considered. (Default = 10 bp)")
	parser.add_option("-r", action="store", dest="search_range", default=5000, type=int, 
						help="Range to search for potentially relevant reads (Default = 5000 bp)")
	parser.add_option("-x", action="store", dest="autosome", default=None, type=int, 
						help="Processes only this autosome. (Default = all are processed)")
	(options, args) = parser.parse_args()
	
	if (len(args)!=3):
		parser.print_help()
		return 1

	# list of autosomes to be processed 
	list_of_autosomes = range(1,23) # default: all autosomes
	if options.autosome is not None:
		list_of_autosomes = [options.autosome]
	
	vcf_filename 		= os.path.abspath(args[0])
	bam_healthy_filename 	= os.path.abspath(args[1])
	bam_cancer_filename 	= os.path.abspath(args[2])

	vcf_reader 			= vcf.Reader(open(vcf_filename))
	bam_healthy_processor 		= BAMProcessor(bam_healthy_filename, search_range = options.search_range, primary_alignments_only = options.primary_only) # processes the BAM file of the control/healthy sample
 	bam_cancer_processor 		= BAMProcessor(bam_cancer_filename, search_range = options.search_range, primary_alignments_only = options.primary_only) # processes the BAM file of the cancer sample

	# Walk through all records in the VCF file
	for vcf_record in vcf_reader:  
		if not returnAutosome(vcf_record.CHROM) in list_of_autosomes: # not the autosome we are currently interested in.
			continue
		if len(vcf_record.ALT) != 1: # records with several alternatives are ignored.
			continue 
		is_deletion 	= isDeletion(vcf_record) 
		is_insertion 	= isInsertion(vcf_record)
		if not (is_deletion or is_insertion): # not an indel
			continue
		if options.deletions_only and is_insertion: # it's a deletion and we consider only insertions
			continue
		if options.insertions_only and is_deletion: # it's an insertion and we consider only deletions
			continue
		if returnIndelLength(vcf_record) < options.min_length: # below the minimal length of an indel.
			continue 

		# Obtain the evidence (internal segment based and overlapping alignments) 
		if is_deletion:
			deletion = Deletion(vcf_record)
			print('-', '\t', returnAutosome(vcf_record.CHROM), '\t', vcf_record.POS, '\t', deletion.length)
			isize, isize_prob, splits, splits_prob = bam_healthy_processor.processDeletion(deletion)
			print('\t'.join(map(str,isize)))
			print('\t'.join(map(str,isize_prob)))
			print('\t'.join(map(str,splits)))
			print('\t'.join(map(str,splits_prob)))
			isize, isize_prob, splits, splits_prob = bam_cancer_processor.processDeletion(deletion)
			print('\t'.join(map(str,isize)))
			print('\t'.join(map(str,isize_prob)))
			print('\t'.join(map(str,splits)))
			print('\t'.join(map(str,splits_prob)))
		else:
			insertion = Insertion(vcf_record)
			print('+', '\t', returnAutosome(vcf_record.CHROM), '\t', vcf_record.POS, '\t', insertion.length)
			isize, isize_prob, splits, splits_prob = bam_healthy_processor.processInsertion(insertion)
			print('\t'.join(map(str,isize)))
			print('\t'.join(map(str,isize_prob)))
			print('\t'.join(map(str,splits)))
			print('\t'.join(map(str,splits_prob)))
			isize, isize_prob, splits, splits_prob = bam_cancer_processor.processInsertion(insertion)
			print('\t'.join(map(str,isize)))
			print('\t'.join(map(str,isize_prob)))
			print('\t'.join(map(str,splits)))
			print('\t'.join(map(str,splits_prob)))
			 
	bam_healthy_processor.close()
	bam_cancer_processor.close()

if __name__ == '__main__':
	sys.exit(main())
