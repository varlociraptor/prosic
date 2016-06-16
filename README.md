# Postprocessing for Somatic Mutation Calling (PROSIC)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/prosic/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/prosic/badges/downloads.svg)](http://bioconda.github.io/recipes/prosic/README.html)

## Installation

### Via Bioconda

PROSIC will soon be available via [Bioconda](https://bioconda.github.io), a distribution
of bioinformatics software for the conda package manager.
Bioconda can be set up in any Linux environment, even without admin rights.
With Bioconda set up, PROSIC can be installed via

  $ conda install prosic

### Compiling from source

PROSIC uses the following C libraries:

* The GNU scientific library ([GSL](http://www.gnu.org/software/gsl/)), in particular:
	* `gsl_math.h`,
	* `gsl_min.h` and
	* `gsl_errno.h`.

* The [GMP](https://gmplib.org/) library for arbitrary precision arithmetic.

Further, the PROSIC depends on the following Python packages:

* [PySAM](https://code.google.com/p/pysam/) for working with BAM/SAM files and
* [PyVCF](https://github.com/jamescasbon/PyVCF) for working with VCF files.

To install PROSIC from source, issue the following commands:

	$ cmake CMakeLists.txt && make
	$ cp bin/prosic-call $BIN_PATH
	$ python setup.py install --user

Where `$BIN_PATH` should be a folder in your `$PATH`, e.g. `~/.local/bin`.
After that, you will have access to the three commands `prosic-call`,
`prosic-extract-observations` and `prosic-annotate`.


## Author

[Louis Dijkstra](https://github.com/louisdijkstra)
