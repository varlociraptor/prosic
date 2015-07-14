#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
import os
import sys
import vcf
import pysam
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__))[:-3] + 'python')

from Indel import *
from Alignments import *
from DefaultBAMProcessor import *
from LaserBAMProcessor import *
from BWABAMProcessor import *

__author__ = "Louis Dijkstra"

usage = """%prog [options] <aligner> <vcf-file> <bam-healthy> <bam-cancer> 

	<aligner>	the aligner that was used. At the moment
			there are three options: 1) laser, 2) bwa
			or 3) default. 
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

def exceedsThreshold (value, threshold): 
	"""Returns whether value exceeds threshold. When threshold is None, 
	   it returns False"""
	if threshold != None: 
		return value > threshold 
	return False

def main():

	parser = OptionParser(usage=usage)
	parser.add_option("--deletions-only", action="store_true", dest="deletions_only", default=False, 
				  		help="Only deletions are stored and insertions are discarded.")
	parser.add_option("--insertions-only", action="store_true", dest="insertions_only", default=False, 
				  		help="Only insertions are stored and deletions are discarded.")
	parser.add_option("--primary-only", action="store_true", dest="primary_only", default=False,
						help="Only primary alignments are taken into account. (Default = False; also secondary etc. alignments are considered.)")
	parser.add_option("-c", action="store", dest="centerpoints_distance_thres_del", default=50, type=int,
				  		help="Maximal difference in centerpoints between split and deletion for the read to be considered evidence for the deletion. Only applicable when the default aligner is used. (Default = 50bp)")
	parser.add_option("-C", action="store", dest="centerpoints_distance_thres_ins", default=50, type=int,
				  		help="Maximal difference in centerpoints between split and insertion for the read to be considered evidence for the insertion. Only applicable when the default aligner is used. (Default = 50bp)")
	parser.add_option("-i", action="store", dest="del_is_threshold", default=None, type=int,
						help="Discards any insert size observations for any deletion that exceeds the given length. (Overwrites default dictated by the aligner used)")
	parser.add_option("-I", action="store", dest="ins_is_threshold", default=None, type=int,
						help="Discards any insert size observations for any insertion that exceeds the given length. (Overwrites default dictated by the aligner used)")
	parser.add_option("-l", action="store", dest="length_thres_del", default=20, type=int,
				  		help="Maximal difference in length between split and deletion for the read to be considered evidence for the deletion. Only applicable when the default aligner is used. (Default = 20bp)")
	parser.add_option("-L", action="store", dest="length_thres_ins", default=20, type=int,
				  		help="Maximal difference in length between split and insertion for the read to be considered evidence for the insertion. Only applicable when the default aligner is used. (Default = 20bp)")
	parser.add_option("-m", action="store", dest="min_length", default=0, type=int,
				  		help="Minimal length of an indel to be considered. (Default = 0 bp, i.e., all)")
	parser.add_option("-r", action="store", dest="search_range", default=5000, type=int, 
						help="Range to search for potentially relevant reads (Default = 5000 bp)")
	parser.add_option("-s", action="store", dest="del_split_threshold", default=None, type=int,
						help="Discards any split read observations for any deletion that exceeds the given length. (Overwrites default dictated by the aligner used)")
	parser.add_option("-S", action="store", dest="ins_split_threshold", default=None, type=int,
						help="Discards any split read observations for any insertion that exceeds the given length. (Overwrites default dictated by the aligner used)")
	(options, args) = parser.parse_args()
	
	if (len(args)!=4):
		parser.print_help()
		return 1
	
	aligner 		= args[0].lower() 
	vcf_filename 		= os.path.abspath(args[1])
	bam_healthy_filename 	= os.path.abspath(args[2])
	bam_cancer_filename 	= os.path.abspath(args[3])

	vcf_reader 		= vcf.Reader(open(vcf_filename))

	# allocate memory
	bam_healthy_processor, bam_cancer_processor = None, None
	if aligner == 'laser': ### LASER ###
		# set defaults (if not set already)
		if options.del_split_threshold == None: options.del_split_threshold 	= 900 
		if options.ins_split_threshold == None: options.ins_split_threshold 	= 60  
		if options.ins_is_threshold == None: 	options.ins_is_threshold 	= 150 		

		bam_healthy_processor = LaserBAMProcessor(bam_healthy_filename, search_range = options.search_range, primary_alignments_only = options.primary_only)
		bam_cancer_processor = LaserBAMProcessor(bam_cancer_filename, search_range = options.search_range, primary_alignments_only = options.primary_only)
	elif aligner == 'bwa' or aligner == 'bwamemm': ### BWA or BWAMEMM ###
		# set defaults (if not set already)
		if options.del_split_threshold == None: options.del_split_threshold 	= 40
		if options.ins_split_threshold == None: options.ins_split_threshold 	= 25  
		if options.ins_is_threshold == None: 	options.ins_is_threshold 	= 150 		

		bam_healthy_processor = BWABAMProcessor(bam_healthy_filename, search_range = options.search_range, primary_alignments_only = options.primary_only)
		bam_cancer_processor = BWABAMProcessor(bam_cancer_filename, search_range = options.search_range, primary_alignments_only = options.primary_only)
	elif aligner == 'default': ### DEFAULT ###
		bam_healthy_processor = DefaultBAMProcessor(bam_healthy_filename, search_range = options.search_range, primary_alignments_only = options.primary_only, centerpoints_thres_del = options.centerpoints_distance_thres_del, centerpoints_thres_ins = options.centerpoints_distance_thres_ins, length_thres_del = options.length_thres_del, length_thres_ins = options.length_thres_ins)
		bam_cancer_processor = DefaultBAMProcessor(bam_cancer_filename, search_range = options.search_range, primary_alignments_only = options.primary_only, centerpoints_thres_del = options.centerpoints_distance_thres_del, centerpoints_thres_ins = options.centerpoints_distance_thres_ins, length_thres_del = options.length_thres_del, length_thres_ins = options.length_thres_ins) 
	else: 
		print('ERROR: aligner %s was not recognized. Options are: laser, bwa or default'%aligner)
		exit() 
		
	# Walk through all records in the VCF file
	for vcf_record in vcf_reader:  
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
			ignore_insert_size_obs 	= exceedsThreshold(deletion.length, options.del_is_threshold)
			ignore_split_obs 	= exceedsThreshold(deletion.length, options.del_split_threshold)
			
			print('-', '\t', vcf_record.CHROM, '\t', vcf_record.POS, '\t', deletion.length)

			# extract observations from the healthy sample 			
			isize, isize_prob, splits, splits_prob = bam_healthy_processor.processDeletion(deletion)
			
			if ignore_insert_size_obs: print('\n') # 2 empty lines   	
			else:				
				print('\t'.join(map(str,isize)))
				print('\t'.join(map(str,isize_prob)))  

			if ignore_split_obs: print('\n') # 2 empty lines
			else: 
				print('\t'.join(map(str,splits)))
				print('\t'.join(map(str,splits_prob)))
	
			# extract observations from the cancer sample 	
			isize, isize_prob, splits, splits_prob = bam_cancer_processor.processDeletion(deletion)
			
			if ignore_insert_size_obs: print('\n') # 2 empty lines   	
			else:				
				print('\t'.join(map(str,isize)))
				print('\t'.join(map(str,isize_prob)))  

			if ignore_split_obs: print('\n') # 2 empty lines
			else: 
				print('\t'.join(map(str,splits)))
				print('\t'.join(map(str,splits_prob)))
		else:
			insertion = Insertion(vcf_record)
			ignore_insert_size_obs 	= exceedsThreshold(insertion.length, options.ins_is_threshold)
			ignore_split_obs 	= exceedsThreshold(insertion.length, options.ins_split_threshold)

			print('+', '\t', vcf_record.CHROM, '\t', vcf_record.POS, '\t', insertion.length)
			
			# extract observations from the healthy sample 			
			isize, isize_prob, splits, splits_prob = bam_healthy_processor.processInsertion(insertion)
			
			if ignore_insert_size_obs: print('\n') # 2 empty lines   	
			else:				
				print('\t'.join(map(str,isize)))
				print('\t'.join(map(str,isize_prob)))  

			if ignore_split_obs: print('\n') # 2 empty lines
			else: 
				print('\t'.join(map(str,splits)))
				print('\t'.join(map(str,splits_prob)))

			# extract observations from the cancer sample 
			isize, isize_prob, splits, splits_prob = bam_cancer_processor.processInsertion(insertion)
			
			if ignore_insert_size_obs: print('\n') # 2 empty lines   	
			else:				
				print('\t'.join(map(str,isize)))
				print('\t'.join(map(str,isize_prob)))  

			if ignore_split_obs: print('\n') # 2 empty lines
			else: 
				print('\t'.join(map(str,splits)))
				print('\t'.join(map(str,splits_prob)))
			 
	bam_healthy_processor.close()
	bam_cancer_processor.close()

if __name__ == '__main__':
	sys.exit(main())
