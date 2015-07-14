----------------------------------------
DEPENDENCIES
----------------------------------------
The application presented here requires the following Python packages to be installed:  

pyvcf - for working with VCF files, see https://github.com/jamescasbon/PyVCF
pysam - for working with BAM/SAM files, see https://github.com/pysam-developers/pysam

----------------------------------------
USAGE
----------------------------------------
The main script 'call_somatic_mutations.py' takes a VCF file with potential somatic variants and two 
BAM files, one for the healthy/control sample and one for the cancer/disease sample, as input. It calls 
every deletion larger than a predefined length (default is 10 bp, see option '-m')  as either 

	(1) somatic - it only appears among the cancer cells;
	(2) germline - it appears among the healthy cells;
	(3) not present - it does not appear among the cancer and healthy cells. 

In addition, it provides a maximum a posteriori (MAP) estimate of the variant allele frequency among the 
healthy and cancer cells. The script outputs an annotated VCF file where the results for every deletion considered
are provided in the INFO field under the identifiers:  

'MAP_CANCER_VAF' 	- the MAP estimate of the cancer cells (varies over the unit interval);
'MAP_HEALTHY_VAF' 	- the MAP estimates of the healthy cells. Value can either be 0.0, 0.5 
			  and 1.0, corresponding to 'absent', 'heterozygous' and 'homozygous'; 
'CALL'			- either 'somatic', 'germline' or 'not present'; the hypthesis with the 
			  highest posterior probability;
'POSTERIOR_PROB' 	- the posterior probability that belong to the call made. 

The script should be used as follows: 

python call_somatic_mutations.py [options] <mu> <sigma> <vcf-file> <bam-healthy> <bam-cancer>  <output-vcf>

where 
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

The script allows for the following optional settings: 

  --align-uncertainty-off
                        Alignment uncertainty is discared. Every alignment is
                        taken for granted.
  --internal-segment-based-only
                        Only internal segment based evidences are taken into
                        account.
  --primary-only        Only primary alignments are taken into account.
                        (Default = False; also secondary etc. alignments are
                        considered.)
  --split-read-alignments-only
                        Only split-read evidences are taken into account.
  --typing-uncertainty-off
                        Typing uncertainty for the overlapping alignments is
                        turned off. Equivalent to setting epsilon_a and
                        epsilon_p both to 0.0, i.e., '-e 0.0 -E 0.0'.
  -a ALPHA              Level of impurity in disease/cancer sample. (Default =
                        0.0; no healthy cells present in the disease/cancer
                        sample)
  -e EPSILON_A          Value of the parameter epsilon_a (the probability of
                        observing a split when the indel is not present;
                        Default = 0.0)
  -E EPSILON_P          Value of the parameter epsilon_p (the probability of
                        not observing a split when the indel is present;
                        Default = 0.0)
  -m MIN_LENGTH         Minimal length of an indel to be considered. (Default
                        = 10 bp)
  -r SEARCH_RANGE       Range to search for potentially relevant reads
                        (Default = 5000 bp)
  -t TOLERANCE          Precision threshold used for numerically maximizing
                        the loglikelihood function and computing integrals.
                        (Default = 10^-8)
  -x AUTOSOME           Processes only this autosome. (Default = all are
                        processed)

----------------------------------------
ASSUMPTIONS
----------------------------------------
The script is implemented under the assumptions that
	- the prior for the healthy and cancer variant allele frequences are uniform;
	- healthy cells are diploid. 

----------------------------------------
CONTACT
----------------------------------------
In case of any questions, feel free to contact me: 

Louis Dijkstra
dijkstra@cwi.nl 


