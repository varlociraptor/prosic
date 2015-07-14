#ifndef SMC_LIKELIHOOD_H_
#define SMC_LIKELIHOOD_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_errno.h>

#include <gmp.h>

#include "smc_input.h" // contains data structure 

#define CONTINUE_MINIMIZATION 	0
#define DONE_MINIMIZATION	1
#define IMPOSSIBLE_HEALTHY_VAF	2 

/* parameters - stay constant throughout all computations */
double ALPHA ; 
double MU_H ; 
double SIGMA_H ; 
double RATE_H ; 
double MU_C ; 
double SIGMA_C ; 
double RATE_C ; 
double EPS_A ; 
double EPS_P ; 
size_t MAX_ITER ; 
double EPSABS ; 
size_t N_PANELS ; 

/* Probability distribution for insert size observations */
double f(double x, double mu, double std, size_t delta, double rate) ; 

/* The PMF of the insert size distribution when there is no indel (null distribution) */
double f0(double x, double mu, double std) ; 

/* Determines the likelihood of h_vaf and c_vaf given the data D */ 
void likelihood (mpf_t* l_final, double h_vaf, double c_vaf, data D) ;

/* Approximates the integral of the likelihood function using the trapezoidal rule */
void integrate_likelihood(mpf_t* I_final, double h_vaf, data D) ; 

/* Approximates the posterior probabilities of the hypotheses 'somatic', 'germline' and 'not present' */
void determinePosteriorProbabilities(double* p_somatic, double* p_germline, double* p_not_present, data D) ; 

/* Loglikelihood function */
double loglikelihood(double h_vaf, double c_vaf, data D) ; 

/* Loglikelihood function where the healthy VAF (h_vaf) is set to 0.0 */
double logl_h_vaf0(double c_vaf, data D) ;

/* Loglikelihood function where the healthy VAF (h_vaf) is set to 0.5 */
double logl_h_vaf1(double c_vaf, data D) ; 

/* Loglikelihood function where the healthy VAF (h_vaf) is set to 1.0 */
double logl_h_vaf2(double c_vaf, data D) ; 

/* Returns 1 when there are very unlikely insert size observations, 0 otherwise. */
size_t unlikelyInsertSize (data D) ;

/* Removes any unlikely insert size observations from the data. */
data removeUnlikelyInsertSizes (data D) ; 

/* Returns 1 when unique global maximum exists, otherwise 0 */
size_t uniqueGlobalMaximumExists (data D) ;

/* Returns the maximum likelihood estimates of h_vaf and c_vaf given data D */
double computeMLE (double* MLE_h_vaf, double* MLE_c_vaf, data D)  ; 

#endif
