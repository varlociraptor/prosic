#ifndef SMC_INPUT_H_
#define SMC_INPUT_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <unistd.h>
#include <float.h>

#define INITIAL_N_OBSERVATIONS 1000
#define LENGTH_CHROM_BUFFER 32 		/* # of characters allocated for chromosome buffer */
#define DELIM	"\t"			/* input file is tab-delimited */

#define VALID_VARIANT 		0 
#define END_OF_FILE_REACHED 	1
#define WRONG_VARIANT_TYPE 	2
#define VARIANT_TOO_SHORT 	3

typedef struct {
	short align_uncertainty_off ; 
	short deletions_only ; 
	short insertions_only ;
	short isize_only ; 
	short split_only ; 

	short verbose ; 

	size_t min_length ; 

	/* Model parameters */
	double alpha ; 		// level of impurity
	double eps_a ; 		// epsilon_a (same for deletions and insertions)
	double eps_p_del ; 	// epsilon_p for deletions
	double eps_p_ins ; 	// epsilon_p for insertions
	double mu_h ; 		// null mean insert sizes healthy sample 
	double sigma_h ;	// null STD insert sizes healthy sample 
	double rate_h ; 	// error rate for the insert size model healthy sample
	double mu_c ; 		// null mean insert sizes cancer sample 
	double sigma_c ; 	// null STD insert sizes cancer sample 
	double rate_c ; 	// error rate for the insert size model cancer sample

	/* parameters for the numerical methods */
	size_t max_iter ; 
	double epsabs ; 
	size_t n_panels ; 
} parameters ; 

/* 
 * Represents the info on the variant
 */
typedef struct {
	char type ; 		// '+' in case of an insertion, '-' in case of a deletion
	char* chromosome ; 	// chromosome that harbours the variant
	size_t position ; 	// position of the variant in the VCF file
	size_t length ; 	// length of the indel 	
} variant ; 

/*
 * Represents the data for one VCF record 
 */
typedef struct {
	// insert size data from the healthy sample
	size_t* h_isize  ; 
	double* h_isize_prob ; 
	size_t h_isize_n ; 
	// split read data from the healthy sample
	size_t* h_split ; 
	double* h_split_prob ; 
	size_t h_split_n ; 
	// insert size data from the cancer sample
	size_t* c_isize ; 
	double* c_isize_prob ; 	
	size_t c_isize_n ; 
	// split read data from the cancer sample
	size_t* c_split ; 
	double* c_split_prob ; 
	size_t c_split_n  ; 
	size_t delta ; 
} data ;

/* Prints usage */
void usage(const char *pname) ; 

/* Parses the options from command line and returns a struct containing all parameters */
void parse_arguments(parameters* p, char* filename, int argc, char **argv) ;

/* Prints the parameters to command line */
void print_parameters(parameters* p) ; 

/* Reads in one line of unknown length */
char *inputString(FILE* fp) ; 

/* Converts a string into an array of integers of unknown length */
size_t *returnIntegerArray(char* buffer, size_t* n_elem) ; 

/* Converts a string into an array of doubles of unknown length */
double *returnDoubleArray(char* buffer, size_t* n_elem) ; 

/* Converts a string into an array of doubles of known length n */
double *returnProbabilities(char* buffer, int n) ; 

/* Returns a list of n doubles initialized as ones */
double *returnListOfOnes(size_t n) ; 

/* Skips the data entries for one variant */
void skipDataEntry (FILE *fp) ; 

/* Returns variant info (when available) */
size_t obtainVariant(FILE *fp, variant *v, parameters *p) ; 

/* Returns the data for one VCF record */
data obtainData(FILE *fp, parameters *p) ;

/* Prints the data to the command line */
void printData (data D) ;  

#endif
