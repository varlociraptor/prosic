/*
** call_somatic_mutations.c
**
**
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "sm_likelihood.h"

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

	size_t i ;
	variant v ; 
	data D ; 

	parameters p = parseOptions(argc, argv) ; 

	while(1) {
		// read in the data 
		if (obtainVariant(fp, &v) == 1) { // eof of file has been reached
			break ; 
		} 
		D = obtainData(fp) ; 

		if (v.type == '+') {
			D.delta = -1.0*v.length ; 
		} else {
			D.delta = v.length ; 
		}
		
		// print info on the VCF record
		printf("%c\t%zd\t%zd\t%zd\n", v.type, v.autosome, v.position, v.length);
		
		
		
		printf("logl: %f\t%f\n", loglikelihood(0.0, 0.0, D), logl_h_vaf0(0.0, D)) ; 
		printf("logl: %f\t%f\n", loglikelihood(0.5, 0.0, D), logl_h_vaf1(0.0, D)) ; 
		printf("logl: %f\t%f\n", loglikelihood(1.0, 0.0, D), logl_h_vaf2(0.0, D)) ; 
	}

    	fclose(fp);
    	exit(EXIT_SUCCESS) ; 
}
