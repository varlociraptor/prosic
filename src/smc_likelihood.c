#include "smc_likelihood.h"

double f(double x, double mu, double std, size_t delta, double rate) {
	/* The PMF of the insert size distribution */
	double d = (double) delta ; 
	double prob = 0.5 / ( rate*(1.0 - 0.5*erfc((mu + 0.5)/std*M_SQRT1_2)) + (1.0 - rate)*(1.0 - 0.5*erfc((mu + d + 0.5)/std* M_SQRT1_2)) ) ;  // normalization factor
	prob *= rate*( erfc((-x - 0.5 + mu)/std*M_SQRT1_2) - erfc((-x + 0.5 + mu)/std*M_SQRT1_2) ) + (1.0 - rate)*( erfc((-x - 0.5 + mu + d)/std*M_SQRT1_2) - erfc((-x + 0.5 + mu + d)/std*M_SQRT1_2) ) ;   
	return prob ; 
}

double f0(double x, double mu, double std) {
	/* The PMF of the insert size distribution when there is no indel (null distribution) */
	double prob = 0.5 / (1.0 - 0.5*erfc((mu + 0.5)/std* M_SQRT1_2)) ;  // normalization factor
	prob *= (erfc((-x - 0.5 + mu)/std*M_SQRT1_2) - erfc((-x + 0.5 + mu)/std*M_SQRT1_2))  ;   
	return prob ; 
}

void likelihood (mpf_t* l_final, double h_vaf, double c_vaf, data D) {
	/* Determines the likelihood of h_vaf and c_vaf given the data D*/
	mpf_t l, obs_l ; // likelihood and likelihood for one observation 
	size_t i ; 

	// initialize and set likelihood
	mpf_init(l) ; 
	mpf_init(obs_l) ; 
	mpf_set_ui(l, 1) ;
	
	for (i = 0; i < D.h_isize_n; i ++) {
		mpf_set_d(obs_l, D.h_isize_prob[i] * (h_vaf*f(D.h_isize[i], MU_H, SIGMA_H, D.delta, RATE_H) + (1.0 - h_vaf)*f0(D.h_isize[i], MU_H, SIGMA_H)) + (1.0 - D.h_isize_prob[i])) ; 
		mpf_mul(l, l, obs_l) ; 
	}

	for (i = 0; i < D.h_split_n; i ++) {
		mpf_set_d(obs_l, D.h_split_prob[i] * (h_vaf*(D.h_split[i] * (1.0 - EPS_P) + (1 - D.h_split[i]) * EPS_P) + (1.0 - h_vaf)*(D.h_split[i]*EPS_A + (1.0 - D.h_split[i])*(1.0 - EPS_A))) + (1.0 - D.h_split_prob[i])) ; 
		mpf_mul(l, l, obs_l) ; 
	}

	for (i = 0; i < D.c_isize_n; i ++) {
		mpf_set_d(obs_l, D.c_isize_prob[i] * (ALPHA * (h_vaf*f(D.c_isize[i], MU_H, SIGMA_H, D.delta, RATE_H) + (1.0 - h_vaf)*f0(D.c_isize[i], MU_H, SIGMA_H)) + (1.0 - ALPHA) * (c_vaf*f(D.c_isize[i], MU_C, SIGMA_C, D.delta, RATE_C) + (1.0 - c_vaf)*f0(D.c_isize[i], MU_C, SIGMA_C))) + (1.0 - D.c_isize_prob[i])) ; 
		mpf_mul(l, l, obs_l) ; 
	}

	for (i = 0; i < D.c_split_n; i ++) {
		mpf_set_d(obs_l, D.c_split_prob[i] * (ALPHA * (h_vaf*(D.c_split[i] * (1.0 - EPS_P) + (1 - D.c_split[i]) * EPS_P) + (1.0 - h_vaf)*(D.c_split[i]*EPS_A + (1.0 - D.c_split[i])*(1.0 - EPS_A))) + (1.0 - ALPHA) * (c_vaf*(D.c_split[i] * (1.0 - EPS_P) + (1 - D.c_split[i]) * EPS_P) + (1.0 - c_vaf)*(D.c_split[i]*EPS_A + (1.0 - D.c_split[i])*(1.0 - EPS_A)))) + (1.0 - D.c_split_prob[i])) ; 
		mpf_mul(l, l, obs_l) ; 
	}

	mpf_set((*l_final), l) ; 
}


void integrate_likelihood(mpf_t* I_final, double h_vaf, data D) {
	/* Approximates the integral of the likelihood function for fixed h_vaf while c_vaf varies over the unit interval. 
	 * We use the trapezoidal rule with N equally spaced panels. 
	 */
	mpf_t I, l, step_size, two; 
	size_t i ; 

	// initialize and set area and likelihood
	mpf_init(I) ; // approximate area 
	mpf_init(l) ; // likelihood
  	mpf_init(step_size) ; 
	mpf_init(two) ; 
	
	mpf_set_d(I, 0) ;
	mpf_set_d(step_size, 1.0 / (double)N_PANELS) ; 
	mpf_set_ui(two, 2) ; 

	likelihood(&l, h_vaf, 0.0, D) ; // compute likelihood at c_vaf = 0	
	mpf_add(I, I, l) ; 
	likelihood(&l, h_vaf, 1.0, D) ; // compute likelihood at c_vaf = 1	
	mpf_add(I, I, l) ; 
	
	// all the points between 0 and 1 
	for (i = 1; i < N_PANELS; i ++) {
		likelihood(&l, h_vaf, i / (double)N_PANELS, D) ;
		mpf_mul(l, l, two) ; // multiply the likelihood with 2
		mpf_add(I, I, l) ; 
	}

	// determine normalization constant (stored in step_size)
	mpf_div(step_size, step_size, two) ; 
	mpf_mul(I, I, step_size) ; // normalize

	mpf_set((*I_final), I) ; 
}

/* Approximates the posterior probabilities of the hypotheses 'somatic', 'germline' and 'not present' */
void determinePosteriorProbabilities(double* p_somatic, double* p_germline, double* p_not_present, data D) {
	mpf_t ps, pg, pnp ; // posterior probabilities (s - somatic, g - germline, np - not present)
	mpf_t I ; // integral

	mpf_init(ps) ; 
	mpf_init(pg) ; 
	mpf_init(pnp) ;
	mpf_init(I) ;

	likelihood (&pnp, 0.0, 0.0, D) ; 
	mpf_div_ui(pnp, pnp, 3) ; 

	integrate_likelihood(&ps, 0.0, D) ; 
	mpf_div_ui(ps, ps, 9) ; 

	integrate_likelihood(&pg, 0.5, D) ;
	integrate_likelihood(&I, 1.0, D) ;
	mpf_add(pg, pg, I) ; 
	mpf_div_ui(pg, pg, 9) ; 

	// compute normalization constant
	mpf_add(I, ps, pg) ; 
	mpf_add(I, I, pnp) ; 

	// normalize probabilities
	mpf_div(ps, ps, I) ; 
	mpf_div(pg, pg, I) ; 
	mpf_div(pnp, pnp, I) ; 

	(*p_somatic) = mpf_get_d(ps) ;   
	(*p_germline) = mpf_get_d(pg) ;
	(*p_not_present) = mpf_get_d(pnp) ;
}


double loglikelihood(double h_vaf, double c_vaf, data D) {
	double logl = 0.0 ;
	size_t i ; 
	// walk through all insert size observations from the healthy sample 
	for (i = 0; i < D.h_isize_n; i ++) {
		logl += log(D.h_isize_prob[i] * (h_vaf*f(D.h_isize[i], MU_H, SIGMA_H, D.delta, RATE_H) + (1.0 - h_vaf)*f0(D.h_isize[i], MU_H, SIGMA_H)) + (1.0 - D.h_isize_prob[i])) ; 
	}
	
	// walk through all split observations from the healthy sample 
	for (i = 0; i < D.h_split_n; i ++) {
		logl += log(D.h_split_prob[i] * (h_vaf*(D.h_split[i] * (1.0 - EPS_P) + (1 - D.h_split[i]) * EPS_P) + (1.0 - h_vaf)*(D.h_split[i]*EPS_A + (1.0 - D.h_split[i])*(1.0 - EPS_A))) + (1.0 - D.h_split_prob[i])) ; 
	}

	// walk through all insert size observations from the cancer sample 
	for (i = 0; i < D.c_isize_n; i ++) {
		logl += log(D.c_isize_prob[i] * (ALPHA * (h_vaf*f(D.c_isize[i], MU_H, SIGMA_H, D.delta, RATE_H) + (1.0 - h_vaf)*f0(D.c_isize[i], MU_H, SIGMA_H)) + (1.0 - ALPHA) * (c_vaf*f(D.c_isize[i], MU_C, SIGMA_C, D.delta, RATE_C) + (1.0 - c_vaf)*f0(D.c_isize[i], MU_C, SIGMA_C))) + (1.0 - D.c_isize_prob[i])) ; 
	}

	// walk through all split observations from the cancer sample
	for (i = 0; i < D.c_split_n; i ++) {
		logl += log(D.c_split_prob[i] * (ALPHA * (h_vaf*(D.c_split[i] * (1.0 - EPS_P) + (1 - D.c_split[i]) * EPS_P) + (1.0 - h_vaf)*(D.c_split[i]*EPS_A + (1.0 - D.c_split[i])*(1.0 - EPS_A))) + (1.0 - ALPHA) * (c_vaf*(D.c_split[i] * (1.0 - EPS_P) + (1 - D.c_split[i]) * EPS_P) + (1.0 - c_vaf)*(D.c_split[i]*EPS_A + (1.0 - D.c_split[i])*(1.0 - EPS_A)))) + (1.0 - D.c_split_prob[i])) ; 
	}

	return logl ; 
}

double logl_h_vaf1(double c_vaf, data D) {
	double logl = 0.0 ;

	size_t i ; 
	// walk through all insert size observations from the healthy sample 
	for (i = 0; i < D.h_isize_n; i ++) {
		logl += log(D.h_isize_prob[i] * (.5*f(D.h_isize[i], MU_H, SIGMA_H, D.delta, RATE_H) + .5*f0(D.h_isize[i], MU_H, SIGMA_H)) + (1.0 - D.h_isize_prob[i])) ; 
	}
	
	// walk through all split observations from the healthy sample 
	for (i = 0; i < D.h_split_n; i ++) {
		logl += log(D.h_split_prob[i] * (.5*(D.h_split[i] * (1.0 - EPS_P) + (1 - D.h_split[i]) * EPS_P) + .5*(D.h_split[i]*EPS_A + (1.0 - D.h_split[i])*(1.0 - EPS_A))) + (1.0 - D.h_split_prob[i])) ; 
	}

	// walk through all insert size observations from the cancer sample 
	for (i = 0; i < D.c_isize_n; i ++) {
		logl += log(D.c_isize_prob[i] * (ALPHA * (.5*f(D.c_isize[i], MU_H, SIGMA_H, D.delta, RATE_H) + .5*f0(D.c_isize[i], MU_H, SIGMA_H)) + (1.0 - ALPHA) * (c_vaf*f(D.c_isize[i], MU_C, SIGMA_C, D.delta, RATE_C) + (1.0 - c_vaf)*f0(D.c_isize[i], MU_C, SIGMA_C))) + (1.0 - D.c_isize_prob[i])) ; 
	}

	// walk through all split observations from the cancer sample
	for (i = 0; i < D.c_split_n; i ++) {
		logl += log(D.c_split_prob[i] * (ALPHA * (.5*(D.c_split[i] * (1.0 - EPS_P) + (1 - D.c_split[i]) * EPS_P) + .5*(D.c_split[i]*EPS_A + (1.0 - D.c_split[i])*(1.0 - EPS_A))) + (1.0 - ALPHA) * (c_vaf*(D.c_split[i] * (1.0 - EPS_P) + (1 - D.c_split[i]) * EPS_P) + (1.0 - c_vaf)*(D.c_split[i]*EPS_A + (1.0 - D.c_split[i])*(1.0 - EPS_A)))) + (1.0 - D.c_split_prob[i])) ; 
	}

	return logl ; 
}


double logl_h_vaf0 (double c_vaf, data D) {
	double logl = 0.0 ; 
	size_t i; 
	// walk through all insert size observations from the healthy sample 
	for (i = 0; i < D.h_isize_n; i ++) {
		logl += log(D.h_isize_prob[i] * f0(D.h_isize[i], MU_H, SIGMA_H) + (1.0 - D.h_isize_prob[i])) ; 
	}
	
	// walk through all split observations from the healthy sample 
	for (i = 0; i < D.h_split_n; i ++) {
		logl += log(D.h_split_prob[i] * (D.h_split[i]*EPS_A + (1.0 - D.h_split[i])*(1.0 - EPS_A)) + (1.0 - D.h_split_prob[i])) ; 
	}

	// walk through all insert size observations from the cancer sample 
	for (i = 0; i < D.c_isize_n; i ++) {
		logl += log(D.c_isize_prob[i] * (ALPHA*f0(D.c_isize[i], MU_H, SIGMA_H) + (1.0 - ALPHA) * (c_vaf*f(D.c_isize[i], MU_C, SIGMA_C, D.delta, RATE_C) + (1.0 - c_vaf)*f0(D.c_isize[i], MU_C, SIGMA_C))) + (1.0 - D.c_isize_prob[i])) ; 
	}

	// walk through all split observations from the cancer sample 
	for (i = 0; i < D.c_split_n; i ++) {
		logl += log(D.c_split_prob[i] * (ALPHA *(D.c_split[i]*EPS_A + (1.0 - D.c_split[i])*(1.0 - EPS_A)) + (1.0 - ALPHA) * (c_vaf*(D.c_split[i] * (1.0 - EPS_P) + (1 - D.c_split[i]) * EPS_P) + (1.0 - c_vaf)*(D.c_split[i]*EPS_A + (1.0 - D.c_split[i])*(1.0 - EPS_A)))) + (1.0 - D.c_split_prob[i])) ; 
	}
	return logl ; 
}

double logl_h_vaf2(double c_vaf, data D) {
	double logl = 0.0 ;
	size_t i ; 
	// walk through all insert size observations from the healthy sample 
	for (i = 0; i < D.h_isize_n; i ++) {
		logl += log(D.h_isize_prob[i] * f(D.h_isize[i], MU_H, SIGMA_H, D.delta, RATE_H) + (1.0 - D.h_isize_prob[i])) ; 
	}
	
	// walk through all split observations from the healthy sample 
	for (i = 0; i < D.h_split_n; i ++) {
		logl += log(D.h_split_prob[i] * (D.h_split[i] * (1.0 - EPS_P) + (1 - D.h_split[i]) * EPS_P) + (1.0 - D.h_split_prob[i])) ; 
	}

	// walk through all insert size observations from the cancer sample 
	for (i = 0; i < D.c_isize_n; i ++) {
		logl += log(D.c_isize_prob[i] * (ALPHA * f(D.c_isize[i], MU_H, SIGMA_H, D.delta, RATE_H) + (1.0 - ALPHA) * (c_vaf*f(D.c_isize[i], MU_C, SIGMA_C, D.delta, RATE_C) + (1.0 - c_vaf)*f0(D.c_isize[i], MU_C, SIGMA_C))) + (1.0 - D.c_isize_prob[i])) ; 
	}

	// walk through all split observations from the cancer sample
	for (i = 0; i < D.c_split_n; i ++) {
		logl += log(D.c_split_prob[i] * (ALPHA * (D.c_split[i] * (1.0 - EPS_P) + (1 - D.c_split[i]) * EPS_P) + (1.0 - ALPHA) * (c_vaf*(D.c_split[i] * (1.0 - EPS_P) + (1 - D.c_split[i]) * EPS_P) + (1.0 - c_vaf)*(D.c_split[i]*EPS_A + (1.0 - D.c_split[i])*(1.0 - EPS_A)))) + (1.0 - D.c_split_prob[i])) ; 
	}

	return logl ; 
}

data removeUnlikelyInsertSizes (data D) {
	/* Remove unlikely insert size observations */
	data newD ; 

	newD.h_split 		= D.h_split ; 
	newD.h_split_prob 	= D.h_split_prob ; 
	newD.h_split_n 		= D.h_split_n ; 
	newD.c_split 		= D.c_split ; 
	newD.c_split_prob 	= D.c_split_prob ; 
	newD.c_split_n 		= D.c_split_n ; 
	newD.delta 		= D.delta ; 

	newD.h_isize_n = 0 ; 
	newD.c_isize_n = 0 ; 

	newD.h_isize 		= malloc(D.h_isize_n * sizeof(size_t)) ; 	
	newD.c_isize 		= malloc(D.c_isize_n * sizeof(size_t)) ; 	
	newD.h_isize_prob 	= malloc(D.h_isize_n * sizeof(double)) ; 	
	newD.c_isize_prob 	= malloc(D.c_isize_n * sizeof(double)) ; 	

	if (newD.h_isize == NULL || newD.c_isize == NULL || newD.h_isize_prob == NULL || newD.c_isize_prob == NULL) {
		printf("ERROR: insufficient memory for allocation.\n") ; 
		exit(EXIT_FAILURE) ; 
	}

	size_t i ; 
	for (i = 0 ; i < D.h_isize_n ; i ++) {
		if (f(D.h_isize[i], MU_H, SIGMA_H, D.delta, RATE_H) != 0.0 || f0(D.h_isize[i], MU_H, SIGMA_H) != 0.0) {
			newD.h_isize[newD.h_isize_n] = D.h_isize[i] ; 
			newD.h_isize_prob[newD.h_isize_n] = D.h_isize_prob[i] ; 
			newD.h_isize_n ++ ; 
		}
	}

	for (i = 0 ; i < D.c_isize_n ; i ++) {
		if (f(D.c_isize[i], MU_C, SIGMA_C, D.delta, RATE_C) != 0.0 || f0(D.c_isize[i], MU_C, SIGMA_C) != 0.0) {
			newD.c_isize[newD.c_isize_n] = D.c_isize[i] ; 
			newD.c_isize_prob[newD.c_isize_n] = D.c_isize_prob[i] ; 
			newD.c_isize_n ++ ; 
		}
	}
	
	return newD ; 
}


size_t unlikelyInsertSize (data D) {
	/* Returns 1 when there are very unlikely insert size observations. 
	 * Unlikely refers here to an insert size for which f returns 0 for both models. 
         */
	size_t i ; 
	for (i = 0; i < D.h_isize_n; i ++) {
		if (f(D.h_isize[i], MU_H, SIGMA_H, D.delta, RATE_H) == 0.0 && f0(D.h_isize[i], MU_H, SIGMA_H) == 0.0) {
			return 1 ; 
		}
	}

	for (i = 0; i < D.c_isize_n; i ++) {
		if (f(D.c_isize[i], MU_C, SIGMA_C, D.delta, RATE_C) == 0.0 && f0(D.c_isize[i], MU_C, SIGMA_C) == 0.0) {
			return 1 ; 
		} 
	}
	return 0 ; 
}
	
size_t uniqueGlobalMaximumExists (data D) {
	/* Returns 1 when unique global maximum exists, otherwise 0 */
	if (ALPHA == 1.0) { // no cancer cells present in the sample
		return 0 ; 
	}

	size_t i ; 
	size_t c_vaf_max_exists = 0 ; // unique global maximum for c_vaf exists (0 = no, 1 = yes)
	if (EPS_A != (1.0 - EPS_P)) {		
		for (i = 0; i < D.c_split_n; i ++) {
			if (D.c_split_prob[i] != 0.0) { 				
				c_vaf_max_exists = 1 ; 
				break ; 
			}
		}
	}		
	for (i = 0; i < D.c_isize_n; i ++) {
		if (D.c_isize_prob[i] != 0.0) {
			if (f0(D.c_isize[i], MU_C, SIGMA_C) != f(D.c_isize[i], MU_C, SIGMA_C, D.delta, RATE_C)) {
				c_vaf_max_exists = 1 ; 
				break ; 
			}
		}
	}

	if (ALPHA != 0.0 || c_vaf_max_exists == 0) {
		return c_vaf_max_exists ; 
	}

	if (EPS_A != (1.0 - EPS_P)) {		
		for (i = 0; i < D.h_split_n; i ++) {
			if (D.h_split_prob[i] != 0.0) { 				
				return 1 ;
			}
		}
	}		
	for (i = 0; i < D.h_isize_n; i ++) {
		if (D.h_isize_prob[i] != 0.0) {
			if (f0(D.h_isize[i], MU_H, SIGMA_H) != f(D.h_isize[i], MU_H, SIGMA_H, D.delta, RATE_H)) {
				return 1 ; 
			}
		}
	}
	return 0 ; 
}

double aux_logl_h_vaf0 (double x, void * params) {
	return -1.0 * logl_h_vaf0(x, (*(data*)params)) ; 
}

double aux_logl_h_vaf1 (double x, void * params) {
	return -1.0 * logl_h_vaf1(x, (*(data*)params)) ; 
}

double aux_logl_h_vaf2 (double x, void * params) {
	return -1.0 * logl_h_vaf2(x, (*(data*)params)) ; 
}

size_t prepareMinimization (double h_vaf, double *x, double *a, double *b, data D) {
	/* Determines whether the h_vaf value is likely at all. In addition, it locates 
	   the unimodal part of the loglikelihood function. */
	mpf_t l0, l1, l2 ; 
	mpf_init(l0) ; 
	mpf_init(l1) ; 
	mpf_init(l2) ; 

	likelihood (&l0, h_vaf, 0.0, D) ; 
	likelihood (&l1, h_vaf, 0.5, D) ; 
	likelihood (&l2, h_vaf, 1.0, D) ; 

	if (mpf_cmp_d(l0, 0.0) == 0 && mpf_cmp_d(l1, 0.0) == 0 && mpf_cmp_d(l2, 0.0) == 0) 
		return IMPOSSIBLE_HEALTHY_VAF ; 

	if (mpf_cmp_d(l0, 0.0) == 0) 
		(*a) = EPSABS / 2 ; 
	else 
		(*a) = 0.0 ; 

	if (mpf_cmp_d(l2, 0.0) == 0) 
		(*b) = 1.0 - EPSABS / 2 ; 
	else 
		(*b) = 1.0 ; 

	/* Find appropriate initial guess x */
	size_t iter = 0 ; 
	double logl_a = loglikelihood(h_vaf, (*a), D); 
	double logl_x = loglikelihood(h_vaf, (*x), D);
	double logl_b = loglikelihood(h_vaf, (*b), D);
	while (iter < MAX_ITER && (logl_x <= logl_a || logl_x <= logl_b)) {
		iter ++ ; 
		if (logl_x <= logl_a) 
			(*x) = ((*x) + (*a)) / 2 ; 
		else 
			(*x) = ((*x) + (*b)) / 2 ; 
		logl_x = loglikelihood(h_vaf, (*x), D) ; 				
	}
	
	if (iter == MAX_ITER) {
		if (logl_a > logl_b) 
			(*x) = (*a) ;
		else 
			(*x) = (*b) ;  
		return DONE_MINIMIZATION ; 
	}
	return CONTINUE_MINIMIZATION ; 
}


double computeMLE (double* MLE_h_vaf, double* MLE_c_vaf, data D) {
	/* Returns the maximum likelihood estimates of h_vaf and c_vaf given data D */
	double mle_c_vaf0 = 0.5, mle_c_vaf1 = 0.5, mle_c_vaf2 = 0.5 ;
	double a, b, max_logl0, max_logl1, max_logl2 ; 
	size_t status, iter = 0 ;
	
	const gsl_min_fminimizer_type *T;
	gsl_min_fminimizer *s;
	T = gsl_min_fminimizer_brent;
  	s = gsl_min_fminimizer_alloc(T);
	gsl_function F;
	F.params = &D ; 

	gsl_set_error_handler_off();		
	
	status = prepareMinimization(0.0, &mle_c_vaf0, &a, &b, D) ; 
	switch(status) {
		case DONE_MINIMIZATION:
			max_logl0 = logl_h_vaf0(mle_c_vaf0, D) ;  break ; 
		case IMPOSSIBLE_HEALTHY_VAF:
			max_logl0 = -INFINITY ; break ; 
		case CONTINUE_MINIMIZATION: 
			iter = 0 ; 
			F.function = &aux_logl_h_vaf0; 
			status = gsl_min_fminimizer_set (s, &F, mle_c_vaf0, a, b);
			do { 
				iter ++ ;
				gsl_min_fminimizer_iterate(s) ;
				status = gsl_min_test_interval(gsl_min_fminimizer_x_lower(s), gsl_min_fminimizer_x_upper(s), EPSABS, 0.0) ; 
			} while (status == GSL_CONTINUE && iter < MAX_ITER);
			mle_c_vaf0 = gsl_min_fminimizer_x_minimum(s) ;
			max_logl0 = -1.0 * gsl_min_fminimizer_f_minimum(s) ; 
			break ; 
	}

	status = prepareMinimization(0.5, &mle_c_vaf1, &a, &b, D) ; 
	switch(status) {
		case DONE_MINIMIZATION:
			max_logl1 = logl_h_vaf1(mle_c_vaf1, D) ; break ; 
		case IMPOSSIBLE_HEALTHY_VAF:
			max_logl1 = -INFINITY ; break ; 
		case CONTINUE_MINIMIZATION: 
			iter = 0 ; 
			F.function = &aux_logl_h_vaf1; 
			status = gsl_min_fminimizer_set (s, &F, mle_c_vaf1, a, b);
			do { 
				iter ++ ;
				gsl_min_fminimizer_iterate(s) ;
				status = gsl_min_test_interval(gsl_min_fminimizer_x_lower(s), gsl_min_fminimizer_x_upper(s), EPSABS, 0.0) ; 
			} while (status == GSL_CONTINUE && iter < MAX_ITER);
			mle_c_vaf1 = gsl_min_fminimizer_x_minimum(s) ;
			max_logl1 = -1.0 * gsl_min_fminimizer_f_minimum(s) ; 
			break ; 
	}

	status = prepareMinimization(1.0, &mle_c_vaf2, &a, &b, D) ; 
	switch(status) {
		case DONE_MINIMIZATION:
			max_logl2 = logl_h_vaf2(mle_c_vaf1, D) ; break ; 
		case IMPOSSIBLE_HEALTHY_VAF:
			max_logl2 = -INFINITY ; break ; 
		case CONTINUE_MINIMIZATION: 
			iter = 0 ; 
			F.function = &aux_logl_h_vaf2; 
			status = gsl_min_fminimizer_set (s, &F, mle_c_vaf2, a, b); 
			do { 
				iter ++ ;
				gsl_min_fminimizer_iterate(s) ;
				status = gsl_min_test_interval(gsl_min_fminimizer_x_lower(s), gsl_min_fminimizer_x_upper(s), EPSABS, 0.0) ; 
			} while (status == GSL_CONTINUE && iter < MAX_ITER);
			mle_c_vaf2 = gsl_min_fminimizer_x_minimum(s) ;
			max_logl2 = -1.0 * gsl_min_fminimizer_f_minimum(s) ; 
			break ; 
	}

	gsl_min_fminimizer_free (s);

	// Determine MLE estimate 
	if (max_logl0 >= max_logl1 && max_logl0 >= max_logl2) {
		(*MLE_h_vaf) = 0.0 ;
		(*MLE_c_vaf) = mle_c_vaf0 ; 
		return max_logl0 ; 
	} 

	if (max_logl1 > max_logl0 && max_logl1 > max_logl2) {
		(*MLE_h_vaf) = 0.5 ;
		(*MLE_c_vaf) = mle_c_vaf1 ; 
		return max_logl1 ; 
	} 

	(*MLE_h_vaf) = 1.0 ;
	(*MLE_c_vaf) = mle_c_vaf2 ; 
	return max_logl2 ; 	 
}


























