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
                
# Het per 10,000 bp

HetDic = defaultdict(list)
HetInput =  str(sys.argv[1])
with gzip.open(HetInput,"rt", newline = '') as file:                                                                                          
        rows = csv.reader(file, delimiter='\t')
        for row in rows:
            if row[0][0] == "#":
                pass
            else:
                PositionBin = int(row[1]) - (int(row[1]) % 100000)
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
HetOutputfile = "HetPer100KB.csv"
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
    

