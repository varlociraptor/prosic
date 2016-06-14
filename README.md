# Postprocessing for Somatic Mutation Calling (PROSIC)


## Installation

### Via Bioconda

PROSIC is available via [Bioconda](https://bioconda.github.io), a distribution
of bioinformatics software for the conda package manager.
Bioconda can be set up in any Linux environment, even without admin rights.
With Bioconda set up, PROSIC can be installed via

  $ conda install prosic

### Compiling from source

PROSIC uses the following C libraries:

* The GNU scientific library (GSL - see http://www.gnu.org/software/gsl/), in particular:
	* `gsl_math.h`,
	* `gsl_min.h` and
	* `gsl_errno.h`.

* The GMP library for arbitrary precision arithmetic, see https://gmplib.org/.

Further, the PROSIC depends on the following Python packages:

* PySAM (see https://code.google.com/p/pysam/) for working with BAM/SAM files and
* PyVCF (see https://github.com/jamescasbon/PyVCF) for working with VCF files.

In order to compile the C-code in the folder `src/`, simply type in the main directory:

	$ cmake CMakeLists.txt && make

The executable `prosic-call` will be placed in the `bin/` folder together with the Python scripts
`prosic-extract-observations` and `prosic-annotate`.


## Author

[Louis Dijkstra](https://github.com/louisdijkstra)
