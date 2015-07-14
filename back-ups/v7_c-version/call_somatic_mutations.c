/*
** call_somatic_mutations.c
**
**
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>

#include "input.h"

#define INITIAL_GUESS 0.5

const char *delim = "\t" ; 

double alpha = 0.0 ; 
double mu_h = 112.0 ; 
double sigma_h = 15.0 ; 
double mu_c = 112.0 ; 
double sigma_c = 15.0 ; 
double eps_a = 0.0 ; 
double eps_p = 0.0 ; 



double f(double x, double mean, double std) {
	return erfc(-((x - mean + 1.0) / std) * M_SQRT1_2) - erfc((mean - x) / std * M_SQRT1_2) ; 
}


double loglikelihood(double h_vaf, double c_vaf, data D) {
	double logl = 0.0 ;
	size_t i = 0 ; 
	// walk through all insert size observations from the healthy sample 
	for (i; i < D.h_isize_n; i ++) {
		logl += log(D.h_isize_prob[i] * (h_vaf*f(D.h_isize[i], mu_h + D.delta, sigma_h) + (1.0 - h_vaf)*f(D.h_isize[i], mu_h, sigma_h)) + (1.0 - D.h_isize_prob[i])) ; 
	}
	
	// walk through all split observations from the healthy sample 
	for (i; i < D.h_split_n; i ++) {
		logl += log(D.h_split_prob[i] * (h_vaf*(D.h_split[i] * (1.0 - eps_p) + (1 - D.h_split[i]) * eps_p) + (1.0 - h_vaf)*(D.h_split[i]*eps_a + (1.0 - D.h_split[i])*(1.0 - eps_a))) + (1.0 - D.h_split_prob[i])) ; 
	}

	// walk through all insert size observations from the cancer sample 
	for (i; i < D.c_isize_n; i ++) {
		logl += log(D.c_isize_prob[i] * (alpha * (h_vaf*f(D.c_isize[i], mu_h + D.delta, sigma_h) + (1.0 - h_vaf)*f(D.c_isize[i], mu_h, sigma_h)) + (1.0 - alpha) * (c_vaf*f(D.c_isize[i], mu_c + D.delta, sigma_c) + (1.0 - c_vaf)*f(D.c_isize[i], mu_c, sigma_c))) + (1.0 - D.c_isize_prob[i])) ; 
	}

	// walk through all split observations from the cancer sample
	for (i; i < D.c_split_n; i ++) {
		logl += log(D.c_split_prob[i] * (alpha * (h_vaf*(D.c_split[i] * (1.0 - eps_p) + (1 - D.c_split[i]) * eps_p) + (1.0 - h_vaf)*(D.c_split[i]*eps_a + (1.0 - D.c_split[i])*(1.0 - eps_a))) + (1.0 - alpha) * (c_vaf*(D.c_split[i] * (1.0 - eps_p) + (1 - D.c_split[i]) * eps_p) + (1.0 - c_vaf)*(D.c_split[i]*eps_a + (1.0 - D.c_split[i])*(1.0 - eps_a)))) + (1.0 - D.c_split_prob[i])) ; 
	}

	return logl ; 
}

double logl_h_vaf1(double c_vaf, data D) {
	double logl = 0.0 ;
	size_t i = 0 ; 
	// walk through all insert size observations from the healthy sample 
	for (i; i < D.h_isize_n; i ++) {
		logl += log(D.h_isize_prob[i] * (.5*f(D.h_isize[i], mu_h + D.delta, sigma_h) + .5*f(D.h_isize[i], mu_h, sigma_h)) + (1.0 - D.h_isize_prob[i])) ; 
	}
	
	// walk through all split observations from the healthy sample 
	for (i; i < D.h_split_n; i ++) {
		logl += log(D.h_split_prob[i] * (.5*(D.h_split[i] * (1.0 - eps_p) + (1 - D.h_split[i]) * eps_p) + .5*(D.h_split[i]*eps_a + (1.0 - D.h_split[i])*(1.0 - eps_a))) + (1.0 - D.h_split_prob[i])) ; 
	}

	// walk through all insert size observations from the cancer sample 
	for (i; i < D.c_isize_n; i ++) {
		logl += log(D.c_isize_prob[i] * (alpha * (.5*f(D.c_isize[i], mu_h + D.delta, sigma_h) + .5*f(D.c_isize[i], mu_h, sigma_h)) + (1.0 - alpha) * (c_vaf*f(D.c_isize[i], mu_c + D.delta, sigma_c) + (1.0 - c_vaf)*f(D.c_isize[i], mu_c, sigma_c))) + (1.0 - D.c_isize_prob[i])) ; 
	}

	// walk through all split observations from the cancer sample
	for (i; i < D.c_split_n; i ++) {
		logl += log(D.c_split_prob[i] * (alpha * (.5*(D.c_split[i] * (1.0 - eps_p) + (1 - D.c_split[i]) * eps_p) + .5*(D.c_split[i]*eps_a + (1.0 - D.c_split[i])*(1.0 - eps_a))) + (1.0 - alpha) * (c_vaf*(D.c_split[i] * (1.0 - eps_p) + (1 - D.c_split[i]) * eps_p) + (1.0 - c_vaf)*(D.c_split[i]*eps_a + (1.0 - D.c_split[i])*(1.0 - eps_a)))) + (1.0 - D.c_split_prob[i])) ; 
	}

	return logl ; 
}


double logl_h_vaf0 (double c_vaf, data D) {
	double logl = 0.0 ; 
	size_t i = 0 ; 
	// walk through all insert size observations from the healthy sample 
	for (i; i < D.h_isize_n; i ++) {
		logl += log(D.h_isize_prob[i] * f(D.h_isize[i], mu_h, sigma_h) + (1.0 - D.h_isize_prob[i])) ; 
	}
	
	// walk through all split observations from the healthy sample 
	for (i; i < D.h_split_n; i ++) {
		logl += log(D.h_split_prob[i] * (D.h_split[i]*eps_a + (1.0 - D.h_split[i])*(1.0 - eps_a)) + (1.0 - D.h_split_prob[i])) ; 
	}

	// walk through all insert size observations from the cancer sample 
	for (i; i < D.c_isize_n; i ++) {
		logl += log(D.c_isize_prob[i] * (alpha*f(D.c_isize[i], mu_h, sigma_h) + (1.0 - alpha) * (c_vaf*f(D.c_isize[i], mu_c + D.delta, sigma_c) + (1.0 - c_vaf)*f(D.c_isize[i], mu_c, sigma_c))) + (1.0 - D.c_isize_prob[i])) ; 
	}

	// walk through all split observations from the cancer sample 
	for (i; i < D.c_split_n; i ++) {
		logl += log(D.c_split_prob[i] * (alpha *(D.c_split[i]*eps_a + (1.0 - D.c_split[i])*(1.0 - eps_a)) + (1.0 - alpha) * (c_vaf*(D.c_split[i] * (1.0 - eps_p) + (1 - D.c_split[i]) * eps_p) + (1.0 - c_vaf)*(D.c_split[i]*eps_a + (1.0 - D.c_split[i])*(1.0 - eps_a)))) + (1.0 - D.c_split_prob[i])) ; 
	}
	return logl ; 
}

double logl_h_vaf2(double c_vaf, data D) {
	double logl = 0.0 ;
	size_t i = 0 ; 
	// walk through all insert size observations from the healthy sample 
	for (i; i < D.h_isize_n; i ++) {
		logl += log(D.h_isize_prob[i] * f(D.h_isize[i], mu_h + D.delta, sigma_h) + (1.0 - D.h_isize_prob[i])) ; 
	}
	
	// walk through all split observations from the healthy sample 
	for (i; i < D.h_split_n; i ++) {
		logl += log(D.h_split_prob[i] * (D.h_split[i] * (1.0 - eps_p) + (1 - D.h_split[i]) * eps_p) + (1.0 - D.h_split_prob[i])) ; 
	}

	// walk through all insert size observations from the cancer sample 
	for (i; i < D.c_isize_n; i ++) {
		logl += log(D.c_isize_prob[i] * (alpha * f(D.c_isize[i], mu_h + D.delta, sigma_h) + (1.0 - alpha) * (c_vaf*f(D.c_isize[i], mu_c + D.delta, sigma_c) + (1.0 - c_vaf)*f(D.c_isize[i], mu_c, sigma_c))) + (1.0 - D.c_isize_prob[i])) ; 
	}

	// walk through all split observations from the cancer sample
	for (i; i < D.c_split_n; i ++) {
		logl += log(D.c_split_prob[i] * (alpha * (D.c_split[i] * (1.0 - eps_p) + (1 - D.c_split[i]) * eps_p) + (1.0 - alpha) * (c_vaf*(D.c_split[i] * (1.0 - eps_p) + (1 - D.c_split[i]) * eps_p) + (1.0 - c_vaf)*(D.c_split[i]*eps_a + (1.0 - D.c_split[i])*(1.0 - eps_a)))) + (1.0 - D.c_split_prob[i])) ; 
	}

	return logl ; 
}


const gsl_min_fminimizer_type *T ; 
T = gsl_min_fminimizer_brent ; 
gsl_min_fminimizer *s ; 
s = gsl_min_fminimizer_alloc(T) ; 
gsl_function F0 ; 
F0.function = &logl_h_vaf0 ; 
const size_t max_iter ; 


double returnMLE_c_vaf0(data D, double* logl) {
	int status ; 	
	int iter = 0 ; 
	F0.params = D ; 
	gsl_min_fminimizer_set(s, &F0, INITIAL_GUESS, 0.0, 1.0) ; 
	
}


int main(int argc, char *argv[])
{
	if(argc != 2) {
		printf("usage: %s <observations-file>\n", argv[0]) ;  
		exit(EXIT_SUCCESS) ; 
	}

	FILE *fp = fopen(argv[1], "r") ; 
	if (fp == NULL) {
		printf("Could not open file %s", argv[1]) ; 
		exit(EXIT_FAILURE) ; 
	}

	char type ; 		// '+' in case of an insertion, '-' in case of a deletion
	size_t autosome ; 	// autosome that harbours the variant
	size_t position ; 	// position of the variant in the VCF file
	size_t length ; 	// length of the indel 
	size_t delta ;

	/* number of observations */
	size_t h_isize_n, h_split_n, c_isize_n, c_split_n ; 

	/* observations */
	size_t *h_isize = malloc(INITIAL_N_OBSERVATIONS * sizeof(size_t)) ; 
	size_t *h_split = malloc(INITIAL_N_OBSERVATIONS * sizeof(size_t)) ; 
	size_t *c_isize = malloc(INITIAL_N_OBSERVATIONS * sizeof(size_t)) ; 
	size_t *c_split = malloc(INITIAL_N_OBSERVATIONS * sizeof(size_t)) ; 
	double *h_isize_prob = malloc(INITIAL_N_OBSERVATIONS * sizeof(double)) ; 
	double *h_split_prob = malloc(INITIAL_N_OBSERVATIONS * sizeof(double)) ; 
	double *c_isize_prob = malloc(INITIAL_N_OBSERVATIONS * sizeof(double)) ; 
	double *c_split_prob = malloc(INITIAL_N_OBSERVATIONS * sizeof(double)) ; 

	if (h_isize == NULL || h_split == NULL || c_isize == NULL || c_split == NULL || h_isize_prob == NULL || h_split_prob == NULL || c_isize_prob == NULL || c_split_prob == NULL) {
		printf("Insufficient memory to allocate the buffer.\n") ;
		exit(EXIT_FAILURE) ; 
	}

	char *buffer ; 
	size_t i ;
	char ch ; 

	while(1) { // read till the eof 
		if(fscanf(fp, "%c%zd%zd%zd%c", &type, &autosome, &position, &length, &ch) != 5){
			break ; 
		}
		

		// print info on the VCF record
		printf("%c\t%zd\t%zd\t%zd\n", type, autosome, position, length);
		
		// collect data 
		data D ; 
		
		if (type == '+') {
			D.delta = -1 * length ; 
		} else {
			D.delta = length ; 
		}

		buffer = inputString(fp);
		D.h_isize = returnIntegerArray(buffer, &D.h_isize_n) ; 
		
		buffer = inputString(fp) ; 
		D.h_isize_prob = returnProbabilities(buffer, D.h_isize_n) ; 

		buffer = inputString(fp);
		D.h_split = returnIntegerArray(buffer, &D.h_split_n) ; 
		
		buffer = inputString(fp) ; 
		D.h_split_prob = returnProbabilities(buffer, D.h_split_n) ; 

		buffer = inputString(fp);
		D.c_isize = returnIntegerArray(buffer, &D.c_isize_n) ; 
		
		buffer = inputString(fp) ; 
		D.c_isize_prob = returnProbabilities(buffer, D.c_isize_n) ; 

		buffer = inputString(fp);
		D.c_split = returnIntegerArray(buffer, &D.c_split_n) ; 
		
		buffer = inputString(fp) ; 
		D.c_split_prob = returnProbabilities(buffer, D.c_split_n) ; 

		printf("logl: %f\t%f\n", loglikelihood(0.0, 0.0, D), logl_h_vaf0(0.0, D)) ; 
		printf("logl: %f\t%f\n", loglikelihood(0.5, 0.0, D), logl_h_vaf1(0.0, D)) ; 
		printf("logl: %f\t%f\n", loglikelihood(1.0, 0.0, D), logl_h_vaf2(0.0, D)) ; 

		//data D = {h_isize, h_isize_prob, h_isize_n, h_split, h_split_prob, h_split_n, c_isize, c_isize_prob, c_isize_n, c_split, c_split_prob, c_split_n, delta} ; 
		//printf("%zd\t%zd\t%zd\n", D.h_isize_n, D.h_split_n, D.c_isize_n) ; 

		


	}


    	fclose(fp);

	// TODO free allocated memory 

    	exit(EXIT_SUCCESS) ; 
}
