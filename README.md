Postprocessing for Somatic Mutation Calling 
===========================================

Dependencies 
------------

The entire project uses the following C libraries: 

* The GNU scientific library (GSL - see http://www.gnu.org/software/gsl/), in particular:
	* `gsl_math.h`,
	* `gsl_min.h` and 
	* `gsl_errno.h`. 

* The GMP library for arbitrary precision arithmetic, see https://gmplib.org/. 

The project depends for Python on the following packages: 

* PySAM (see https://code.google.com/p/pysam/) for working with BAM/SAM files and
* PyVCF (see https://github.com/jamescasbon/PyVCF) for working with VCF files. 
* snakemake (see https://bitbucket.org/johanneskoester/snakemake/wiki/Home).  

Installation
------------

In order to compile the C-code in the folder `src/`, simply type in the main directory: 

	$ cmake CMakeLists.txt && make

The executable `sm_caller` will be placed in the `bin/` folder together with the Python scripts.   

Internal Reports
----------------

The repository contains the directory `internal-reports/`. This folder contains a number of reports on the methods involved in the project. Here follows a list of its subdirectories and a short description of the report/info it contains: 

* `confidence-intervals/` - describes how one can obtain confidence intervals (using the likelihood ratio method).

* `integrating-likelihood-function/` - derivation of the integral of the posterior distribution. An O(N^2) algorithm is provided for computing the coefficients of the resulting n+1 order polynomial. 
 
* `non-nested-models/` - discusses whether the models are non-nested (somatic/germline or absent)

* `preprint/` - contains the preprint for the RECOMB paper. Contains the proof that a unique global maximum exists for the likelihood function (!).

* `sm-parameter-space/` - visualization of the two-dimensional parameter space for the latent variable model.

* `somatic-mutation-calling-methods/` - short description of the methods involved. (Currently out-of-date)
