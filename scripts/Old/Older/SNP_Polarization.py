#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 14:49:02 2023

@author: goreill1
"""

import csv
import gzip
import sys
import re
import statistics
import os
import csv
import numpy as np
from collections import defaultdict

#ElephantPCA_NoChendra.afreq
#Elephant_Ancestral_SNPs.afreq

AncestalSNP = str(snakemake.input[0])
NonAncestralSNP = str(snakemake.input[1])

with open(AncestalSNP, 'r') as file:
    line_reader = csv.reader(file, delimiter='\t')
    for line in line_reader:

        if line[0] != "#CHROM":
            if line[0].replace(".", "_") in locals (): # IF THE CHROMDICTONARY EXISTS
                #ADD TO THE DICTONARY IN THE FORMAT: CHROM= {GENEID: TYPE, START, END, LENGTH}
                
                eval(line[0].replace(".", "_"))[line[1].split(":")[1]] += [line[3]]
    
            else: # IF THE CHROMDICTONARY DOSEN'T EXISTS
                exec(line[0].replace(".", "_") + "= defaultdict(list)") #MAKE A DICTONARY NAMED AFTER THE CHROMOSOME 
                #ADD TO THE DICTONARY IN THE FORMAT: CHROM= {GENEID: TYPE, START, END, LENGTH}
                eval(line[0].replace(".", "_"))[line[1].split(":")[1]] += [line[3]]



Population_SNPs = []
with open(NonAncestralSNP, 'r') as file:
    line_reader = csv.reader(file, delimiter='\t')
    for line in line_reader:
        if line[0] == "#CHROM":
            Population_SNPs.append(line+["Polarization"])
        elif line[0][0:2] == "NW":
            pass
        else:
            if line[0].replace(".", "_") in locals ():
                if line[1].split(":")[1] in eval(line[0].replace(".", "_")):
                    Population_SNPs.append([line[0],line[1],line[2],line[3],str(1-float(line[4])),line[5],"Polarized"])
                    print("ZZZZap!")
                else:
                    Population_SNPs.append([line[0],line[1],line[2],line[3],line[4],line[5],"Not_polarized"])
            else:
                pass


Outputfile = str(snakemake.output)
with open(Outputfile, 'w') as file:
    file.writelines('\t'.join(i) + '\n' for i in Population_SNPs)
    