import csv
import gzip
import sys
import re
import statistics
import os
import csv
import numpy as np
from collections import defaultdict


VCF_file = str(sys.argv[1])
polar = -1
OutputData = [["CHROM","POS","REF","ALT","ALT_FREQS","Samples","Homo1_Count","Homo0_Count","Het_Count","Pol"]]

with gzip.open(VCF_file, 'rt') as file:
    line_reader = csv.reader(file, delimiter='\t')
    for line in line_reader:

        Homo1 = 0
        Homo0 = 0
        Hetero = 0
        Polarized = 0

        if line[0][0] == "#":
            if line[0] == "#CHROM":
                for i in range(len(line)):
                    if line[i].find('ERR2260497') >= 0:
                        Polar = i
                        print(i)

        elif len(line[3]) > 1:
            pass
        else:
            if str(line[Polar][0:3]) == '1/1' or str(line[Polar][0:3]) == '1|1' or str(line[Polar][0:3]) == '0/1' or str(line[Polar][0:3]) == '1/0' or str(line[Polar][0:3]) == '0|1' or str(line[Polar][0:3]) == '1|0':
                Polarized = 1

            for geno in line[9:]:
                if geno != line[Polar]:
                    if str(geno[0:3]) == '1/1' or str(geno[0:3]) == '1|1':
                        Homo1 += 1 
                    elif str(geno[0:3]) == '0/0' or str(geno[0:3]) == '0|0':
                        Homo0 += 1
                    elif str(geno[0:3]) == '0/1' or str(geno[0:3]) == '0|1' or str(geno[0:3]) == '1/0' or str(geno[0:3]) == '1|0' :
                        Hetero += 1
            AltCount = (Hetero+Homo1)

            #header is "CHROM","POS","REF","ALT","ALT_FREQS","Samples","Homo1_Count","Homo0_Count","Het_Count"
            OutputData.append([line[0],line[1],line[3],line[4],str(AltCount),str(len(line)-9),str(Homo1),str(Homo0),str(Hetero),str(Polarized)])


Output_file = str(sys.argv[2])

with open(Output_file, 'w') as file:
    for i in OutputData:
        file.writelines('\t'.join(i) + '\n')
