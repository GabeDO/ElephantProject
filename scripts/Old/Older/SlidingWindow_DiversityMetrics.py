#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 15:51:11 2023

@author: goreill1
"""

import gzip
import sys
import re
import statistics
import os
import csv
import numpy as np
import json
import math
from collections import defaultdict


# ues: /Users/goreill1/Desktop/Elephant_SNPs_3000Sample.vcf


                
def GetAlleleProportions(vcf_LociDATA_list):
    p1_count = 0
    p2_count = 0
    Missing = 0
    Homozygotes_count = 0
    Heterozygotes_count = 0
    
    for d in vcf_LociDATA_list:
        if d[0] == '.' or d[2] == '.':
            Missing += 1
        elif d[0] == d[2] and d[0] != '.':
            Homozygotes_count += 1
        else:
            Heterozygotes_count +=1
         
        if d[0] == '0':
            p1_count += 1 
        elif d[2] == '0':
            p1_count += 1
            
        if d[0] == '1':
            p2_count += 1 
        elif d[2] == '1':
            p2_count += 1
        
    if p1_count+p2_count != 0:
        
        p2_freq = p2_count/(p1_count+p2_count)
        
        p1_freq = p1_count/(p1_count+p2_count)
        
        Het_Obs = Heterozygotes_count/(Heterozygotes_count+Homozygotes_count)
        D_0, D_1, D_2 = SNP_QProfile(p1_freq,p2_freq)
        
        return Het_Obs, D_0, D_1, D_2, Missing
    else:
        return 0, 0, 1, 2, Missing


def SNP_QProfile(p1,p2):
    if p1 == 0 or p2 == 0:
        return 1, 1, 1
    else:
        for q in [0,1,2]:
            if q != 1:
                exec("global D" + str(q)+"; D" + str(q) + " = ((p1**q)+(p2**q))**(1/(1-q))" )
            else:
                exec("global D" + str(q)+"; D" + str(q) + " = math.exp( -( (p1*math.log(p1))+(p2*math.log(p2)) ))" )
        
        return D0, D1, D2
    
            
        
            
#data = "/Users/goreill1/Desktop/Elephant_SNPs_NoMissingNoChendra.vcf.gz"
#data = "/Users/goreill1/Desktop/Elephant_SNPs_3000Sample.vcf.gz"
data = sys.argv[1]

with gzip.open(data, "rt", newline = '') as file: 
    rows = csv.reader (file, delimiter='\t')
    for i in rows:
        if i[0][0] == "#" and len(i) > 1:
            NumberOfSamples = len(i)-9
        elif i[0][0] != "#" :
            #i[0] = Chromosome
            #i[1] = position
            #i[5] = quality
            #i[9:] = data
            
            if i[0].replace(".", "_") in locals(): #If chromDicExists
                eval(i[0].replace(".", "_"))[i[1]]  += GetAlleleProportions(i[9:])
            else:
                exec(i[0].replace(".", "_") + "= defaultdict(list)")
                eval(i[0].replace(".", "_"))[i[1]]  += GetAlleleProportions(i[9:])
                
         
                
##iterate through all chromosome dicts
maxLen = 0
OverlapValue = 5000
WindowSize = 10000
CurrentLocal = list(locals())
SlidingWindowData = defaultdict(list)
for i in CurrentLocal:
    if i[0:2] == "NC":
        MaxValue =int(list(eval(i))[-1])
        for x in range(MaxValue):
            if x % WindowSize == 0:
                Het_obs = []
                for key, values in eval(i).items():
                    if int(key) - x > -OverlapValue and x - int(key) < OverlapValue:
                        Het_obs += [eval(i)[key][0]]
                SlidingWindowData[i] += [[x,sum(Het_obs)/len(Het_obs)]]


with open('TEST_SlidingWindow.txt', 'w') as DataSet:
    write = csv.writer(DataSet)
    write.writerow(["Chromosome", "MiddleWindow", "Obs_Het"])
    for key, value in SlidingWindowData.items():
        for v in value:
            write.writerows([[key,v[0],v[1]]])
