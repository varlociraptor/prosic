# Postprocessing for Somatic Mutation Calling (PROSIC)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/prosic/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/prosic/badges/downloads.svg)](http://bioconda.github.io/recipes/prosic/README.html)

## Installation

PROSIC is available via [Bioconda](https://bioconda.github.io), a distribution
of bioinformatics software for the conda package manager.
Bioconda can be set up in any Linux environment, even without admin rights.
With Bioconda set up, PROSIC can be installed via

	$ conda install prosic

## Usage

PROSIC offers three command line utilities:

* `prosic-extract-observations`
* `prosic-call`
* `prosic-annotate`

The purpose of PROSIC is to call somatic insertions and deletions (indels) on tumor/normal sample pairs.
For this, PROSIC requires a VCF file with preliminary indel calls, e.g. obtained with [Platypus](http://www.well.ox.ac.uk/platypus).
Then, calling with PROSIC consists of three steps.

### Step 1: extracting observations

First, we extract informations about relevant variants and the corresponding alignments from
the VCF and the corresponding BAM files:

	$ prosic-extract-observations $ALIGNER pre-calls.vcf tumor.bam normal.bam > observations.txt

See `prosic-extract-observations --help` for information about the valid values of
`$ALIGNER` and other parameters.

### Step 2: calling

Second, the significance of mutations is assessed via

	$ prosic-call observations.txt > calls.txt

See `prosic-call -h` for information about parameters.

### Step 3: annotation

Finally, the initial VCF file is annotated with the somatic mutation calls from PROSIC:

	$ prosic-annotate pre-calls.vcf calls.txt > calls.vcf

Thereby, each indel in the resulting VCF is classified into `CALL=SOMATIC` or `CALL=GERMLINE` and the
posterior probability for the selected event is provided in the field `POSTERIOR_PROB`.

## Compiling from source

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
