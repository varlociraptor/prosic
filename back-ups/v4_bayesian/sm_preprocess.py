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

%prog [options] --truth <bam-healthy> <bam-cancer>

Otherwise: 

%prog [options] <vcf-file> <bam-healthy> <bam-cancer> 

Preprocesses the given VCF file and two BAM files (control and disease/cancer) for calling somatic mutations. 

	<vcf-file> 	tabix-indexed VCF file containing the 
				indels which allele frequency 
				needs to be estimated. 
				NOTE: File needs to be sorted 
				and indexed. 
	<bam-healthy>	BAM file containing the 
				alignments of the healthy sample. 
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
Every line represents an indel of at least a certain minimal length (see option '-m'). The data points are seperated by tabs. 

The first six columns provide info on the deletion: 

type			- 'DEL' for a deletion and 'INS' for an insertion	
chromsome 		- the chromosome on which the indel(is thought to) occur(s)
position   		- reference position as given in the VCF file 
length			- the length of the indel
true healthy vaf	- the true variant allele frequency of healthy cells (can either be 0.0, 0.5 or 1.0)
true cancer vaf		- the true variant allele frequency of the cancer cells 

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

# 24 empty lists for the 22 autosomes and the two sex-chromosomes X and Y  
list_true_deletions = [[] for i in range(24)] 
list_true_insertions = [[] for i in range(24)] 
   
def returnIndex (chromosome):
	if len(chromosome) > 3 and chromosome[:3] == 'chr':
		chromosome = chromosome[3:]
	if chromosome == 'X' or chromosome == 'x':
		return 22
	if chromosome == 'Y' or chromosome == 'y':
		return 23	
	return int(chromosome)

def obtainListOfTrueIndels (min_length, only_deletions, only_insertions): 
	"""Fills the global lists list_true_deletions and list_true_insertions with the true indels from the true VCF file."""
	for true_vcf_record in VCF_TRUTH_READER:
		is_deletion = isDeletion(true_vcf_record)
		is_insertion = isInsertion(true_vcf_record)
		if not is_deletion and not is_insertion: # not an indel
			continue 	
		if returnIndelLength(true_vcf_record) < min_length: # too short
			continue
		if not only_insertions and is_deletion:
			list_true_deletions[returnIndex(true_vcf_record.CHROM)].append(true_vcf_record)
		if not only_deletions and is_insertion:
			list_true_insertions[returnIndex(true_vcf_record.CHROM)].append(true_vcf_record)


def retrieveMatchingTrueDeletion (vcf_record, distance_threshold, length_threshold):
	"""Returns a matching vcf_record representing a deletion from the true VCF file. If there is no similar one, 'None' is returned."""
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

def retrieveMatchingTrueInsertion (vcf_record, distance_threshold, length_threshold):
	"""Returns a matching vcf_record representing an insertion from the true VCF file. If there is no similar one, 'None' is returned."""
	index = returnIndex(vcf_record.CHROM)
	insertion = Insertion(vcf_record)
	candidates = [] # list with potentially matching candidates
	for true_vcf_record in list_true_insertions[index]:
		if abs(true_vcf_record.POS - vcf_record.POS) > 250:
			continue
		true_insertion = Insertion(true_vcf_record)
		if abs(insertion.position - true_insertion.position) <= distance_threshold and abs(insertion.length - true_insertion.length) <= length_threshold:
			return true_vcf_record
	return None
		
def determineVAF (gt_nums):
	if gt_nums == '1|1': 	return 1.0 
	elif gt_nums == None: 	return 0.0 
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
	   for determining the indels (then it is easy). When false, we first need to locate the corresponding truely present genetic variant and 
	   retrieve the truth. The coverage determines with what frequency certain cancer populations are present."""
	if true_vcf_file_used: 
		return returnTrueVAFs(coverage, vcf_record)
	else: # other VCF than the truth is used
		true_vcf_record = None # allocate memory
		if isDeletion(vcf_record):
			true_vcf_record = retrieveMatchingTrueDeletion (vcf_record, distance_threshold, length_threshold) 
		else:
			true_vcf_record = retrieveMatchingTrueInsertion (vcf_record, distance_threshold, length_threshold) 

		if true_vcf_record is None:
			return 0.0, 0.0 # indel does not appear 
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

def main():

	parser = OptionParser(usage=usage)
	parser.add_option("--low", action="store_true", dest="low_precision", default=False,
						help="Uses low precision when matching detected indels with the indels in the true VCF file, i.e., the centerpoints must be less or exactly 100 bp away and the length must not differ more than 100 bp. Tihs option overwrites the options '-d' and '-l' (!).")
	parser.add_option("--deletions-only", action="store_true", dest="deletions_only", default=False,
						help="Only deletions are processed.")
	parser.add_option("--insertions-only", action="store_true", dest="insertions_only", default=False,
						help="Only insertions are processed.")
	parser.add_option("--primary-only", action="store_true", dest="primary_only", default=False,
						help="Only primary alignments are taken into account. (Default = False)")
	parser.add_option("--truth", action="store_true", dest="truth", default=False,
						help="The indels in the VCF file containing the true genetic variants are called (/data1/structvar/vaf-experiments/vcf/cancerclones.vcf.gz). No VCF file should be provided. This file is normally only used for validating the results on another VCF file.")
	parser.add_option("-d", action="store", dest="distance_threshold", default=10, type=int,
				  		help="Distance between centerpoints of deletions threshold. The distance must be smaller than this in order for the deletions to be deemed similar (in addition, the lengths must be similar, see option '-l'). (Default = 10 bp)")
	parser.add_option("-l", action="store", dest="length_threshold", default=20, type=int,
				  		help="Lengths of deletions threshold. The length difference must be smaller than this in order for the deletions to be deemed similar (in addition, the placement of the centerpoints must be similar, see option '-d'). (Default = 20 bp)")
	parser.add_option("-m", action="store", dest="min_length", default=10, type=int,
				  		help="Minimal length of an indel to be considered. (Default = 10 bp)")
	parser.add_option("-r", action="store", dest="search_range", default=5000, type=int, 
						help="Range to search for potentially relevant reads (Default = 5000 bp)")
	(options, args) = parser.parse_args()
	
	vcf_filename, bam_healthy_filename, bam_cancer_filename = None, None, None # allocate memory

	if options.truth:
		if (len(args)!=2):
			parser.print_help()
			return 1
		vcf_filename 		= VCF_TRUTH_FILENAME
		bam_healthy_filename 	= os.path.abspath(args[0])
		bam_cancer_filename 	= os.path.abspath(args[1])
	else:
		if (len(args)!=3):
			parser.print_help()
			return 1
		vcf_filename 		= os.path.abspath(args[0])
		bam_healthy_filename 	= os.path.abspath(args[1])
		bam_cancer_filename 	= os.path.abspath(args[2])

		# determine the thresholds used for determining whether an indel is similar to another
		if options.low_precision:
			options.length_threshold, options.distance_threshold = 100, 100

		obtainListOfTrueIndels(options.min_length, options.deletions_only, options.insertions_only)	

	coverage = determineCoverage(bam_cancer_filename) # can either be 40 or 80. Is needed for determining the ground truth since it differs from both files
	
	vcf_reader 			= vcf.Reader(open(vcf_filename))
	bam_healthy_processor 		= BAMProcessor(bam_healthy_filename, search_range = options.search_range, primary_alignments_only = options.primary_only) # processes the BAM file of the control/healthy sample
 	bam_cancer_processor 		= BAMProcessor(bam_cancer_filename, search_range = options.search_range, primary_alignments_only = options.primary_only) # processes the BAM file of the cancer sample

	# Walk through all records in the VCF file
	for vcf_record in vcf_reader:  
		if len(vcf_record.ALT) != 1: # records with several alternatives are ignored
			continue 
		is_deletion = isDeletion(vcf_record)
		is_insertion = isInsertion(vcf_record)
		if not is_deletion and not is_insertion: # not an indel
			continue 	
		if returnIndelLength(vcf_record) < options.min_length: # too short
			continue
		if not options.insertions_only and is_deletion: # process the deletion
			print('DEL\t', end = "")			
			deletion = Deletion(vcf_record)
			deletion.print()
			true_h_vaf, true_c_vaf = obtainTruth(vcf_record, options.truth, coverage, options.distance_threshold, options.length_threshold) # obtains the true VAF of the healthy and cancer cells
			print(true_h_vaf, '\t', true_c_vaf, '\t', end = "") # Print the truth 
			# Obtain the evidence (internal segment based and overlapping alignments): 
			h_paired_end_alignments, h_overlapping_alignments = bam_healthy_processor.processDeletion(deletion)
			c_paired_end_alignments, c_overlapping_alignments = bam_cancer_processor.processDeletion(deletion)

			printRawData(h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments) # Prints the raw data and alignment probabilities to the screen	
			print('') # newline
		if not options.deletions_only and is_insertion: # process insertion
			print('INS\t', end = "")			
			insertion = Insertion(vcf_record)
			insertion.print()
			true_h_vaf, true_c_vaf = obtainTruth(vcf_record, options.truth, coverage, options.distance_threshold, options.length_threshold) # obtains the true VAF of the healthy and cancer cells
			print(true_h_vaf, '\t', true_c_vaf, '\t', end = "") # Print the truth 
			# Obtain the evidence (internal segment based and overlapping alignments): 
			h_paired_end_alignments, h_overlapping_alignments = bam_healthy_processor.processInsertion(insertion)
			c_paired_end_alignments, c_overlapping_alignments = bam_cancer_processor.processInsertion(insertion)
			printRawData(h_paired_end_alignments, h_overlapping_alignments, c_paired_end_alignments, c_overlapping_alignments) # Prints the raw data and alignment probabilities to the screen	
			print('') # newline
		
	bam_healthy_processor.close()
	bam_cancer_processor.close()

if __name__ == '__main__':
	sys.exit(main())
