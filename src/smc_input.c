#include "smc_input.h"

void usage(const char *pname) {
	printf(	"Usage: %s [OPTION] <filename.observations>\n"
		"\n"
		" -a\tNUM\tLevel of impurity in the cancer sample (alpha).\n"
		" -e\tNUM\tProbability of a split presence when a deletion is there.\n"
		" -E\tNUM\tProbability of a splits presence when an insertion is there.\n"
		" -y\tNUM\tProbability of a splits absence when an indel is there.\n"
		" -m\tNUM\tMean of insert sizes when the indel is absent (both samples).\n"
		" -s\tNUM\tSTD of insert sizes when the indel is absent (both samples).\n"
		" -m\tNUM\tMean of insert sizes when the indel is absent (both samples).\n"
		" -r\tNUM\tError rate for the insert size model (for both samples).\n"
		" -f\tNUM\tMean of insert sizes when the indel is absent (healthy sample).\n"
		" -g\tNUM\tSTD of insert sizes when the indel is absent (healthy sample).\n"
		" -q\tNUM\tError rate for the insert size model (healthy sample).\n"
		" -F\tNUM\tMean of insert sizes when the indel is absent (cancer sample).\n"
		" -G\tNUM\tSTD of insert sizes when the indel is absent (cancer sample).\n"
		" -Q\tNUM\tError rate for the insert size model (cancer sample).\n"
		" -l\tNUM\tLength of indel to be considered.\n"
		" -P\tNUM\tPrecision used for maximization.\n"
		" -M\tNUM\tMax. number of iterations used for maximization.\n"
		" -N\tNUM\tNumber of panels used for integration (trapezoidal rule).\n"
		" -B\t\tOnly split read evidence.\n"
		" -d\t\tOnly deletions.\n"
		" -i\t\tOnly insertions.\n"
		" -I\t\tOnly insert size evidence.\n"
		" -u\t\tAlignment uncertainty is neglected.\n"
		" -v\t\tVerbose.\n"
		" -h\t\tPrint this help.\n\n", 
		pname) ; 
	exit(0) ; 
}

void print_parameters(parameters* p) {
	printf(	"===PARAMETER SETTINGS===\n\n"
		"-------model parameters-------\n"
		"alpha\t\t%.2f\n"	
		"eps_a\t\t%.2f\n"
		"eps_p_del\t%.2f\n"
		"eps_p_ins\t%.2f\n"
		"mu_h\t\t%.2f\n"
		"sigma_h\t\t%.2f\n"
		"rate_h\t\t%.2f\n"
		"mu_c\t\t%.2f\n"
		"sigma_c\t\t%.2f\n"
		"rate_c\t\t%.2f\n"
		"------------------------------\n\n"
		"-------options-------\n"
		"align. uncertainty off:\t%zd\n" 
		"deletions only:\t\t%zd\n" 
		"insertions only:\t%zd\n" 
		"insert sizes only:\t%zd\n" 
		"split reads only:\t%zd\n"
		"---------------------\n\n"
		"-------numerics-------\n"
		"max_iter\t%zd\n"
		"epsabs\t\t%f\n"
		"n_panels\t%zd\n"
		"----------------------\n\n", 
		p->alpha, p->eps_a, 
		p->eps_p_del, p->eps_p_ins, 
		p->mu_h, p->sigma_h, p->rate_h,
		p->mu_c, p->sigma_c, p->rate_c,
		p->align_uncertainty_off, 
		p->deletions_only, 
		p->insertions_only,
		p->isize_only,
		p->split_only,
		p->max_iter,
		p->epsabs,
		p->n_panels
		) ; 
}

void parse_arguments(parameters* p, char* input_filename, int argc, char **argv) {
	
	/* set defaults */
	p->align_uncertainty_off 	= 0 ; 
	p->deletions_only 		= 0 ; 
	p->insertions_only		= 0 ; 
	p->isize_only			= 0 ; 
	p->split_only			= 0 ; 
	p->verbose			= 0 ; 
	p->min_length			= 0 ;
	p->alpha 			= 0.0 ; 
	p->eps_a 			= 0.0001 ;
	p->eps_p_del			= 0.0961 ; 
	p->eps_p_ins			= 0.3457 ;  
	p->mu_h 			= 112.0 ;
	p->sigma_h			= 15.0 ;
	p->rate_h			= 0.0 ; 
	p->mu_c 			= 112.0 ;
	p->sigma_c			= 15.0 ;
	p->rate_c			= 0.0 ; 
	p->max_iter 			= 100 ; 
	p->epsabs			= 0.0001 ; 
	p->n_panels			= 100 ; 

	size_t ch ; 
	while ((ch = getopt(argc, argv, "a:e:E:y:m:s:r:f:F:g:G:h:H:l:P:M:N:BdiIuvh")) != -1)
    	{
       		switch(ch) {
			case 'a': p->alpha = strtod(optarg, 0); break ; 
			case 'e': p->eps_p_del = strtod(optarg, 0); break ; 
			case 'E': p->eps_p_ins = strtod(optarg, 0); break ; 
			case 'y': p->eps_a = strtod(optarg, 0); break ; 
			case 'm': p->mu_h = strtod(optarg, 0); p->mu_c = p->mu_h; break ; 
			case 's': p->sigma_h = strtod(optarg, 0); p->sigma_c = p->sigma_h; break ; 
			case 'r': p->rate_h = strtod(optarg, 0); p->rate_c = p->rate_h; break ; 
			case 'f': p->mu_h = strtod(optarg, 0); break ; 
			case 'F': p->mu_c = strtod(optarg, 0); break ; 
			case 'g': p->sigma_h = strtod(optarg, 0); break ; 
			case 'G': p->sigma_c = strtod(optarg, 0); break ; 
			case 'q': p->rate_h = strtod(optarg, 0); break ; 
			case 'Q': p->rate_c = strtod(optarg, 0); break ; 
			case 'l': p->min_length = strtol(optarg, 0, 10); break ; 
			case 'P': p->epsabs = strtod(optarg, 0); break ; 
			case 'M': p->max_iter = strtol(optarg, 0, 10); break ; 
			case 'N': p->n_panels = strtol(optarg, 0, 10); break ; 
			case 'B': p->split_only = 1; break ;  
			case 'd': p->deletions_only = 1; break ; 
			case 'i': p->insertions_only = 1; break ;
			case 'I': p->isize_only = 1; break ; 
			case 'u': p->align_uncertainty_off = 1; break ; 
			case 'v': p->verbose = 1; break ; 
        		case 'h': default: usage(argv[0]);
        	}
    	}

	if (argv[optind] == NULL) {
		usage(argv[0]) ; 
		exit(EXIT_SUCCESS) ; 
	}
	
	// check whether parameter settings are valid 

	size_t error_occurred = 0 ; 

	if (p->alpha < 0.0 || p->alpha >= 1.0) {
		printf("ERROR: Invalid argument. Alpha (-a) must lie in [0,1).\n") ;
		error_occurred = 1 ;  
	}

	if (p->eps_p_del < 0.0 || p->eps_p_del >= 1.0) {
		printf("ERROR: Invalid argument. epsilon_p for deletions (-e) must lie in [0,1).\n") ;
		error_occurred = 1 ;  
	}

	if (p->eps_p_ins < 0.0 || p->eps_p_ins >= 1.0) {
		printf("ERROR: Invalid argument. epsilon_p for insertions (-E) must lie in [0,1).\n") ;
		error_occurred = 1 ;  
	}

	if (p->eps_a < 0.0 || p->eps_a >= 1.0) {
		printf("ERROR: Invalid argument. epsilon_a (-y) must lie in [0,1).\n") ;
		error_occurred = 1 ;   
	}

	if (p->mu_h < 0.0) {
		printf("ERROR: Invalid argument. mu_h (-m or -f) must be nonnegative.\n") ;
		error_occurred = 1 ;  
	}

	if (p->mu_c < 0.0) {
		printf("ERROR: Invalid argument. mu_c (-m or -F) must be nonnegative.\n") ;
		error_occurred = 1 ;  
	}

	if (p->sigma_h <= 0.0) {
		printf("ERROR: Invalid argument. sigma_h (-s or -g) must be positive.\n") ;
		error_occurred = 1 ;  
	}

	if (p->sigma_c <= 0.0) {
		printf("ERROR: Invalid argument. sigma_c (-s or -G) must be positive.\n") ;
		error_occurred = 1 ;  
	}

	if (p->rate_h < 0.0 || p->rate_h >= 1.0) {
		printf("ERROR: Invalid argument. rate_h (-r or -q) must lie in [0, 1).\n") ;
		error_occurred = 1 ;  
	}	

	if (p->rate_c < 0.0 || p->rate_c >= 1.0) {
		printf("ERROR: Invalid argument. rate_c (-r or -Q) must lie in [0, 1).\n") ;
		error_occurred = 1 ;  
	}	

	if (p->min_length < 0) {
		printf("ERROR: Invalid argument. min_length (-l) must be nonnegative.\n") ;
		error_occurred = 1 ;  
	}

	if (p->max_iter < 1) {
		printf("ERROR: Invalid argument. max_iter (-M) must be at least 1.\n") ;
		error_occurred = 1 ;  	
	}

	if (p->n_panels < 1) {
		printf("ERROR: Invalid argument. n_panels (-N) must be at least 1.\n") ;
		error_occurred = 1 ;  	
	}

	if (p->epsabs <= 0.0) {
		printf("ERROR: Invalid argument. Precision (-P) must be positive.\n") ;
		error_occurred = 1 ;  
	}

	if (error_occurred == 1) {
		exit(EXIT_FAILURE) ; 
	}

	sprintf(input_filename, argv[optind]) ;
}

char *inputString(FILE* fp) {
	/* Reads in one line of unknown length */
	size_t buffer_length = BUFSIZ ;
	char *buffer = malloc(buffer_length * sizeof(char)) ; 
	if (buffer == NULL) {
		printf("Insufficient memory to allocate the buffer.\n") ;
		exit(EXIT_FAILURE) ; 
	}

	int ch ; 
	size_t i = 0; 
	while (EOF != (ch=fgetc(fp)) && ch != '\n') {
		buffer[i++] = ch ; 
		if (i == buffer_length) {
			buffer = realloc(buffer, sizeof(char) * (buffer_length += BUFSIZ)) ; 
			if (buffer == NULL) {
				printf("Insufficient memory to reallocate the buffer.\n") ;
				exit(EXIT_FAILURE) ; 
			}
		}
	}
	buffer[i++] = '\0' ; 
	return realloc(buffer, sizeof(char) * i); 
}

size_t *returnIntegerArray(char* buffer, size_t* n_elem) { 
	/* Converts a string into an array of integers of unknown length */
	(*n_elem) = 0 ; 
	size_t len = INITIAL_N_OBSERVATIONS ; 
	size_t *observations = malloc(len * sizeof(size_t)) ; 
 	if (observations == NULL) {
		printf("Insufficient memory to allocate the buffer.\n") ;
		exit(EXIT_FAILURE) ; 
	}
	
	size_t i = 0 ; 
	char *token = strtok(buffer, DELIM) ; 
	while (token != NULL) {
		observations[i++] = (size_t) atoi(token) ;
		(*n_elem)++ ; 
		if (i == len) {
			observations = realloc(observations, sizeof(size_t) * (len += INITIAL_N_OBSERVATIONS)) ; 
			if (observations == NULL) {
				printf("Insufficient memory to reallocate the buffer.\n") ;
				exit(EXIT_FAILURE) ; 
			}
		} 
		token = strtok(NULL, DELIM) ; 
		
	}
	observations[i++] = '\0' ; 
	return realloc(observations, sizeof(size_t) * i) ; 	
}

double *returnDoubleArray(char* buffer, size_t* n_elem) { 
	/* Converts a string into an array of doubles of unknown length */
	(*n_elem) = 0 ; 
	size_t len = INITIAL_N_OBSERVATIONS ; 
	size_t *observations = malloc(len * sizeof(double)) ; 
 	if (observations == NULL) {
		printf("Insufficient memory to allocate the buffer.\n") ;
		exit(EXIT_FAILURE) ; 
	}
	
	size_t i = 0 ; 
	
	char *token = strtok(buffer, DELIM) ; 
	while (token != NULL) {
		sscanf(token, "%lf", observations + (i++)) ; 
		(*n_elem)++ ; 
		if (i == len) {
			observations = realloc(observations, sizeof(double) * (len += INITIAL_N_OBSERVATIONS)) ; 
			if (observations == NULL) {
				printf("Insufficient memory to reallocate the buffer.\n") ;
				exit(EXIT_FAILURE) ; 
			}
		} 

		printf("%s\n", token) ; 
		token = strtok(NULL, DELIM) ; 
		
	}
	observations[i++] = '\0' ; 
	return realloc(observations, sizeof(double) * i) ; 	
}

double *returnProbabilities(char* buffer, int n) {
	/* Converts a string into an array of doubles of known length n */
	double *probabilities = malloc(n * sizeof(double)) ; 
	if (probabilities == NULL) {
		printf("Insufficient memory to allocate the buffer.\n") ;
		exit(EXIT_FAILURE) ; 
	}
	
	char *token = strtok(buffer, DELIM) ;
	size_t i = 0 ; 
	for (i; i < n; i ++) {
		sscanf(token, "%lf", probabilities + i) ;
		token = strtok(NULL, DELIM) ; 
	}
	return probabilities ; 
}

double *returnListOfOnes(size_t n) {
	/* Returns a list of n doubles initialized as ones */
	double *list = malloc(n * sizeof(double)) ; 
	if (list == NULL) {
		printf("Insufficient memory to allocate the buffer.\n") ;
		exit(EXIT_FAILURE) ; 
	}
	
	size_t i ; 
	for (i = 0; i < n; i ++) {
		list[i] = 1.0 ; 
	}
	return list ; 
}

size_t obtainVariant(FILE *fp, variant *v, parameters *p) {
	/* Returns variant info and stores it in v 
	 * In case there is no variant anymore, the function returns END_OF_FILE_REACHED.
	 * When the variant is of the wrong type (e.g., only deletions are wanted while
	 * the variant is an insertion), WRONG_VARIANT_TYPE is returned. 
	 * When the variant is too short, it returns VARIANT_TOO_SHORT. 
	 * When the variant is valid, it outputs VALID_VARIANT. 
	 */
	char ch ; 
	if(fscanf(fp, "%c%s%zd%zd%c", &(v->type), v->chromosome, &(v->position), &(v->length), &ch) != 5){
			printf("End of file reached\n") ; 
			return END_OF_FILE_REACHED ; 
	}
	if (p->deletions_only == 1) {
		if (v->type == '+') {return WRONG_VARIANT_TYPE ;} // an insertion, while only deletions are considered
	}
	if (p->insertions_only == 1) {
		if (v->type == '-') {return WRONG_VARIANT_TYPE ;} // a deletion, while only insertions are considered  
	}
	if (p->min_length > v->length) {return VARIANT_TOO_SHORT ;} // too short
	if (v->type == '*') {return WRONG_VARIANT_TYPE ;} // not an indel
	return VALID_VARIANT ;
}

void skipLine (FILE *fp) {
	char c ; 
	do { 
		c = fgetc(fp);
	} while (c != '\n') ; 
}

void skipDataEntry (FILE *fp) {
	size_t i = 0 ;
	while (i < 8) {
		if (fgetc(fp) == '\n') {
			i ++ ; 
		}
	}
}

data obtainData(FILE *fp, parameters *p) {
	/* Returns the data for one VCF record */
	data D ; 
	char *buffer ; 

	// read in the insert sizes from the healthy sample 
	if (p->split_only) {
		D.h_isize_n = 0 ; 
		skipLine(fp) ; 
		skipLine(fp) ; 
	} else {
		buffer = inputString(fp); 
		D.h_isize = returnIntegerArray(buffer, &D.h_isize_n) ; 
		if (p->align_uncertainty_off) {
			skipLine(fp) ; 
			D.h_isize_prob = returnListOfOnes(D.h_isize_n) ; 
		} else {
			buffer = inputString(fp) ;
			D.h_isize_prob = returnProbabilities(buffer, D.h_isize_n) ; 
		}
	}

	// read in the split read observations from the healthy sample 
	if (p->isize_only) {
		D.h_split_n = 0 ; 
		skipLine(fp) ; 
		skipLine(fp) ; 
	} else {
		buffer = inputString(fp); 
		D.h_split = returnIntegerArray(buffer, &D.h_split_n) ; 
		if (p->align_uncertainty_off) {
			skipLine(fp) ; 
			D.h_split_prob = returnListOfOnes(D.h_split_n) ; 
		} else {
			buffer = inputString(fp) ;
			D.h_split_prob = returnProbabilities(buffer, D.h_split_n) ; 
		}
	}


	// read in the insert sizes from the cancer sample 
	if (p->split_only) {
		D.c_isize_n = 0 ; 
		skipLine(fp) ; 
		skipLine(fp) ; 
	} else {
		buffer = inputString(fp); 
		D.c_isize = returnIntegerArray(buffer, &D.c_isize_n) ; 
		if (p->align_uncertainty_off) {
			skipLine(fp) ; 
			D.c_isize_prob = returnListOfOnes(D.c_isize_n) ; 
		} else {
			buffer = inputString(fp) ;
			D.c_isize_prob = returnProbabilities(buffer, D.c_isize_n) ; 
		}
	}

	// read in the split read observations from the cancer sample 
	if (p->isize_only) {
		D.c_split_n = 0 ; 
		skipLine(fp) ; 
		skipLine(fp) ; 
	} else {
		buffer = inputString(fp); 
		D.c_split = returnIntegerArray(buffer, &D.c_split_n) ; 
		if (p->align_uncertainty_off) {
			skipLine(fp) ; 
			D.c_split_prob = returnListOfOnes(D.c_split_n) ; 
		} else {
			buffer = inputString(fp) ;
			D.c_split_prob = returnProbabilities(buffer, D.c_split_n) ; 
		}
	}
	return D;
} 


void printData (data D) {
	/* Prints the data to the command line */
	size_t i ; 
	
	printf("healthy sample\n") ; 
	printf("insert sizes | #obs. %zd\n", D.h_isize_n) ; 
	for (i = 0; i < D.h_isize_n ; i ++) {
		printf("%zd/%.2f ", D.h_isize[i], D.h_isize_prob[i]) ; 
	}
	printf("\n") ;
	printf("split alignm | #obs. %zd\n", D.h_split_n) ; 
	for (i = 0; i < D.h_split_n ; i ++) {
		printf("%zd/%.2f ", D.h_split[i], D.h_split_prob[i]) ; 
	} 
	printf("\n") ;
	printf("cancer sample\n") ; 
	printf("insert sizes | #obs. %zd\n", D.c_isize_n) ; 
	for (i = 0; i < D.c_isize_n ; i ++) {
		printf("%zd/%.2f ", D.c_isize[i], D.c_isize_prob[i]) ; 
	}
	printf("\n") ;
	printf("split alignm | #obs. %zd\n", D.c_split_n) ; 
	for (i = 0; i < D.c_split_n ; i ++) {
		printf("%zd/%.2f ", D.c_split[i], D.c_split_prob[i]) ; 
	} 
	printf("\n") ;
}

