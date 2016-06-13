import pysam
import random

random.seed(46572)

bam = pysam.AlignmentFile(snakemake.input[0], "rb")

subsamples = [pysam.AlignmentFile(f, "wb", template=bam) for f in snakemake.output]

for read in bam:
    random.choice(subsamples).write(read)

for f in subsamples:
    f.close()
