#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 12:32:37 2023

@author: goreill1
"""
import gzip
import sys
import re
import statistics
import os
import csv
import numpy as np
from collections import defaultdict

csv_columns = ['Chromosome','Position']



            
                
# calculate FIS
FisDic = defaultdict(list)
FisList = []
swag = 0

FisInput = str(snakemake.input[0])
with gzip.open(FisInput,"rt", newline = '') as file:                                                                                          
        rows = csv.reader(file, delimiter='\t')
        for row in rows:
            p = 0
            q = 0
            if row[0][0] == "#":
                if row[0] == "#CHROM":
                    for i in row[9:]:
                        ElephantID = i.split("/")
                        csv_columns.append(ElephantID[0])
                else:
                    pass
            elif "INDEL" in  row[7]:
                swag += 1
                pass
            else:
                Loci = str(row[0])+ ", " + str(row[1])+ ", "
                if Loci not in FisDic:
                    FisDic[Loci] = [[],[]]
                
                for Elephant in range(len(row[9:])):
                    #HetObs
                    if row[9+Elephant][0] == "." or row[9+Elephant][2] == ".":
                        pass
                    elif row[9+Elephant][0] == row[9+Elephant][2]:
                        FisDic[Loci][0] += [0]
                    elif row[9+Elephant][0] != row[9+Elephant][2]:
                        FisDic[Loci][0] += [1]
                    
                    #HetExp
                    if row[9+Elephant][0] == "0":
                        p += 1
                    elif row[9+Elephant][0] == "1":
                        q += 1
                    if row[9+Elephant][2] == "0":
                        p += 1
                    elif row[9+Elephant][2] == "1":
                        q += 1
                
                if len(FisDic[Loci][0]) > 0:
                    HetObs = sum(FisDic[Loci][0])/len(FisDic[Loci][0])
                else:
                    pass
                    
                if p == 0 and q ==0:
                    pass
                else:
                    pp = p/(p+q)
                    qq = q/(p+q)
                    HetExp = 1 - (pp**2 + qq**2)
                
                if HetExp == 0:
                    Fis = 0
                else:    
                    Fis = (HetExp - HetObs) / HetExp
                
                FisList.append(Fis)
                
TotalFis = np.mean(FisList)

FisOutputfile = str(snakemake.output[1])
with open(FisOutputfile, 'w') as f:
    f.write('Fis = ' + str(TotalFis))
                    
                
                    
                
# Het per 10,000 bp

HetDic = defaultdict(list)
HetInput = str(snakemake.input[1])
with gzip.open(HetInput,"rt", newline = '') as file:                                                                                          
        rows = csv.reader(file, delimiter='\t')
        for row in rows:
            if row[0][0] == "#":
                pass
            else:
                PositionBin = int(row[1]) - (int(row[1]) % 10000)
                ID = str(row[0])+ ", " + str(PositionBin)+ ", "
                
                if len(HetDic[ID]) == 0:
                    for Elephant in row[9:]:
                        HetDic[ID] += [[]]
                for Elephant in range(len(row[9:])):
                    if row[9+Elephant][0] == "." or row[9+Elephant][2] == ".":
                        pass
                    elif row[9+Elephant][0] == row[9+Elephant][2]:
                        HetDic[ID][Elephant] += [0]
                    elif row[9+Elephant][0] != row[9+Elephant][2]:
                        HetDic[ID][Elephant] += [1]

for key, values in HetDic.items():
    for i in range(len(HetDic[key])):
        if len(HetDic[key][i]) > 0:
            HetDic[key][i] = sum(HetDic[key][i])/len(HetDic[key][i])
        else:
            HetDic[key][i] = "nan"


#Output Het per 10kb
HetOutputfile = str(snakemake.output[0])
csv_file = HetOutputfile 
with open(csv_file, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows([csv_columns])
    
    for key, values in HetDic.items():
        RowOutput = [key.split(",")[0]]
        RowOutput.append(int(key.split(",")[1]))
        for i in values:
            RowOutput.append(i)
        writer.writerows([RowOutput])
    

