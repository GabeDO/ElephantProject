import allel
import scipy
from collections import defaultdict
import gzip
import csv

path_vcf = snakemake.input[0]

CHROM_LIST_NOSEX = snakemake.params[0]

pilist = []

for Chromosome in CHROM_LIST_NOSEX:
    data_vcf = allel.read_vcf(path_vcf, region=Chromosome)
    
    GT = allel.GenotypeArray(data_vcf['calldata/GT'])
    ac = GT.count_alleles()

    
    pi, windows, n_bases, counts = allel.windowed_diversity(data_vcf['variants/POS'], ac, size=10000, step=1000)
    
    for i in range(len(pi)):
        pilist.append([Chromosome, pi[i], windows[i][0], n_bases[i], counts[i]])


OutputFile = str(snakemake.output[0])

with open(OutputFile, "w") as f:
    writer = csv.writer(f)
    writer.writerow(["Chromosome","Pi","Window","Bases","counts"]) 
    writer.writerows(pilist)   
