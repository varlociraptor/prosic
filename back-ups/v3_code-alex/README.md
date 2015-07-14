Post-Processing for Somatic Mutation Calling 
============================================

The code section in our GitHub repository contains now all necessary files. 

The process is divided into two steps to make it easier to experiment. 

Preprocessing the VCF and BAM files
-----------------------------------

A preprocessing step that transforms a given VCF file and two BAM files (one of the healthy and one of the cancer BAM) into a file that contains most relevant information, i.e., the MLEs, confidence intervals and the raw data with their respective alignment probabilities. This preprocessing step is done with

	$ python test_somatic_mutation_caller.py <mu> <sigma> <vcf-file> <healthy-BAM> <cancer-BAM> > output_file.txt

NOTE: look in the help for the extensive list of options! One can set the valeu of alpha (level of impurity), the epsilons, add a cap etc. 

I already generated quite some intermediate results. They can be found in `/data1/structvar/vaf-experiments/vcf/results`

The file `result_files_metadata.ods` contains all relevant information with which settings these files were generated (cap/alpha/epsilons etc.)

All other Python scripts in our code folder works with these intermediate files, e.g., plotting scripts to summarizing the results. 

Generating the results
----------------------

After the preprocessing step, one can generate plots and summarize the results. Here, the final calls (e.g., somatic/germline/not present) are made. 

The main script is `summarize_results.py`, i.e., 

	$ python summarize_results.py r.low.40x.softcip.2-01.txt

It outputs the tables for the genotype results, the main test (calling somatic/germline/not present) and others. Again, there are quite some options, such as the p-value threshold that is used for calling. So make sure you set them appropriately. Particularly, note that the option `--table` provides you with the possibility to examine individual 3x3 tables. 

Plotting the results
--------------------

All files that start with `plot_` generate a plot and can only be run on your home computer (running it on Walter will result in an error message since there is no screen it can send the output to). The next sections we will discuss shortly their use. 

### Creating Genotyping plots

The command 

	$ python plot_genotyping_results.py r.low.40x.softclip.2-01.txt

will generate a genotype plot in the style of the MATE-CLEVER paper. The left column represents the deletions that we call 'not present', the center the ones that we think are heterozygous and the right column are those called homozygous. The colors represent the truth. (With option `-c` one can control the confidence level.)

### Creating calling plots (somatic/germline/not present)

The command 

	$ python plot_calling_results.py r.low.40x.softclip.2-01.txt

will generate a plot as in the previous section, but now for the classes 'Somatic', 'Germline' and 'Not present', i.e., the left column represents the deletions that we call 'somatic', the center the ones that we think are 'germline' and the right column are those called 'not present'. The colors represent the truth. 

NOTE: The options are quite extensive and important to set. For example, 

	$ python plot_calling_results.py -p 0.01 --bonf r.low.40x.softclip.2-01.txt

The p-value treshold is now 1% (in stead of the default of 5%) and a bonferroni correction is applied. 

### Creating maximum likelihood estimate plots  

There are two commands: 

	$ python plot_cvaf_mles.py r.low.40x.softclip.2-01.txt

and 

	$ python plot_vaf_estimates_cancer_bam.py r.low.40x.softclip.2-01.txt

The first one creates a boxplot of the cancer VAF MLEs set against the true VAF of the cancer clones and the cancer clones only. The latter generate a plot of the VAF of the whole cancer BAM file taking both healthy and cancer cells into account. 

Setting alpha (the level of impurity) makes a different for the first command, but not for the second. It could be of interest to look at

	$ python plot_cvaf_mles.py r.low.40x.softclip.2-01.txt

and 

	$ python plot_cvaf_mles.py r.low.40x.a25.softclip.2-01.txt 

(Alpha is 25%, i.e., the truth. `a25`)

It will show that we can estimate the VAF of the cancer clones quite well, when alpha is chosen to be (close to) the truth.

## Important info 

The function that is probably of greatest interest to you is in `ResultsInstance.py`. The class `ResultsInstance` contains the function
`call()` which makes the actual decision; somatic/germline/not present. 







