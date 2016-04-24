import sys
import os

#methods = ['platypus']
#methods = ['manta']
methods = ['pindel']
datasets = ['dream6.dream6']
modes = ['split-alon']
alphas = ['0.0']
#methodnames = ['prosic1']
#methodnames = ['prosic2']
#methodnames = ['prosic3']
methodnames = ['prosic4']
#minmapps = ['0.0','0.9','0.99','0.999']
#minmapps = ['1.0']
minmapps = ['0.0']

chromosomes = ['2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']

for method in methods:
    for dataset in datasets:
        for mode in modes:
            for alpha in alphas:
                for methodname in methodnames:
                    for minmapp in minmapps:

                        command = "cat %s/dreamsubmission/%s.%s.alpha-%s.%s-%s.chr1.vcf " % (method, dataset, mode, alpha, methodname, minmapp)

                        for chromosome in chromosomes:
                            command += "%s/dreamsubmission/%s.%s.alpha-%s.%s-%s.chr%s.noheader.vcf " % (method, dataset, mode, alpha, methodname, minmapp, chromosome)

                        command += "> %s/dreamsubmission/%s.%s.alpha-%s.%s-%s.vcf " % (method, dataset, mode, alpha, methodname, minmapp)
                        
                        os.system(command)
