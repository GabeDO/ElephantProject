#!/usr/bin/env python3 
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 11:58:11 2023

@author: goreill1
"""
 

import gzip
import sys
import re
import statistics
import os
import subprocess
import csv
import numpy as np
import math
##from statsmodels.tsa.stattools import adfuller
from collections import defaultdict

#TandemDups_all.csv
CoverageOutputCounter = 0

SVData_input = str(snakemake.input[0])
SV_TYPE = str(snakemake.params[0])
PolarizationElephant = str(snakemake.params[1])
SVData_Genotyped_Output = str(snakemake.output[0])
sampleloc = str(snakemake.params[3])+"/"


def SetUpCovPerChromForElephant(Elephant_,Chromosome_,path__):
    Process_ = "samtools depth -r " + str(Chromosome_) + " " +  str(path__) + "|  awk '{sum+=$3} END { print sum/NR}' "

    
    if str(Elephant_)+"_CovDICT" in globals():
        Cov = subprocess.check_output(Process_, shell=True)
        eval(str(Elephant_)+"_CovDICT")[Chromosome_] += float(Cov) 
    else:
        exec("globals()['"+ str(Elephant_)+"_CovDICT']" + "= defaultdict(int)") 
        Cov = subprocess.check_output(Process_, shell=True)
        eval(str(Elephant_)+"_CovDICT")[Chromosome_] += float(Cov)

def SetUpChangePointDictForElephant(Elephant_,Chrom ,StartPOS ,EndPOS): #NOT CURRENTLY USED
    SVinCPA = False
    if str(Elephant_.replace('-', ''))+str(Chrom.replace('-', ''))+"_CPADICT" not in locals():
        Process_ = "samtools depth -aa "+sampleloc+str(Elephant_)+str(Elephant_)+"_sampe_sorted.bam | /projects/rogers_research/Code/ChangePointYao"
        ps = subprocess.Popen(Process_,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        psoutput = ps.communicate()[0]
        for line in psoutput:
            ID = str(line[1])
            exec(str(Elephant_.replace('-', ''))+"_CPADICT[ID] = [line[2]/line[3],line[4],line[5]]")
    
    for key, values in eval(str(Elephant_.replace('-', ''))+str(Chrom.replace('-', ''))+"_CPADICT").items():
        if StartPOS >= values[2] and StartPOS <= values[3] and EndPOS >= values[2] and EndPOS <= values[3]:
            SVinCPA = True
            break

    return SVinCPA



ExtraRange = 0
def GetHeterozygosity(data_, read_type):
    PTally = 0
    MTally = 0

    SV_CoverageData = []
    Chromosome = data_[0]

    if read_type == "Rearrangements":
        covindex = 9
        startPOS = min(int(data_[5]),int(data_[6]))-ExtraRange
        endPOS = max(int(data_[5]),int(data_[6]))+ExtraRange
        if startPOS < 1:
            startPOS = 1

        Chromosome_m = data_[2]
        if Chromosome_m == "=":
            Chromosome_m = Chromosome
        startPOS_m = min(int(data_[7]),int(data_[8]))-ExtraRange
        endPOS_m = max(int(data_[7]),int(data_[8]))+ExtraRange
        if startPOS_m < 1:
            startPOS_m = 1


    elif read_type == "TandemDups":
        covindex = 5
        startPOS = min(int(data_[3]),int(data_[4]))-ExtraRange
        endPOS = max(int(data_[3]),int(data_[4]))+ExtraRange
        if startPOS < 1:
            startPOS = 1

    ListOfElephants = ElephantNamesToList(data_[-1])

    for ELEPHANT in ListOfElephants:
        if ELEPHANT != PolarizationElephant and data_[covindex] != "Fail": 
            BAM_PATH =  sampleloc + str(ELEPHANT) + "/" + str(ELEPHANT) + "_sampe_sorted.bam"
            CHROM_REGION = str(str(Chromosome) + ":" + str(startPOS) + "-" + str(endPOS)).replace(" ", "")

            ProcessToRun_SVcov = "samtools depth -r " + CHROM_REGION + " -b NonTE_Regions.bed " +  str(BAM_PATH) + "|  awk '{sum+=$3} END { print sum/NR}' "

            if read_type == "Rearrangements":
                CHROM_REGION_m = str(str(Chromosome_m) + ":" + str(startPOS_m) + "-" + str(endPOS_m)).replace(" ", "")
                ProcessToRun_SVcov_m = "samtools depth -r " + CHROM_REGION_m + " -b NonTE_Regions.bed " +  str(BAM_PATH) + "|  awk '{sum+=$3} END { print sum/NR}' "

            try:
                try:
                    if Chromosome in eval(str(ELEPHANT.replace('-', ''))+"_CovDICT"):
                        TotalCov = eval(str(ELEPHANT.replace('-', ''))+"_CovDICT")[Chromosome]
                    else:
                        SetUpCovPerChromForElephant(ELEPHANT.replace('-', ''),Chromosome,BAM_PATH)
                        TotalCov = eval(str(ELEPHANT.replace('-', ''))+"_CovDICT")[Chromosome]
                except NameError:
                    SetUpCovPerChromForElephant(ELEPHANT.replace('-', ''),Chromosome,BAM_PATH)
                    TotalCov = eval(str(ELEPHANT.replace('-', ''))+"_CovDICT")[Chromosome]
            except NameError:
                SetUpCovPerChromForElephant(ELEPHANT.replace('-', ''),Chromosome,BAM_PATH)
                TotalCov = eval(str(ELEPHANT.replace('-', ''))+"_CovDICT")[Chromosome]
            
            if SV_TYPE == "Rearrangements":
                try:
                    try:
                        if Chromosome_m in eval(str(ELEPHANT.replace('-', ''))+"_CovDICT"):
                            TotalCov_m = eval(str(ELEPHANT.replace('-', ''))+"_CovDICT")[Chromosome_m]
                        else:
                            SetUpCovPerChromForElephant(ELEPHANT.replace('-', ''),Chromosome_m,BAM_PATH)
                            TotalCov_m = eval(str(ELEPHANT.replace('-', ''))+"_CovDICT")[Chromosome_m]
                    except NameError:
                        SetUpCovPerChromForElephant(ELEPHANT.replace('-', ''),Chromosome_m,BAM_PATH)
                        TotalCov_m = eval(str(ELEPHANT.replace('-', ''))+"_CovDICT")[Chromosome_m]
                except NameError:
                    SetUpCovPerChromForElephant(ELEPHANT.replace('-', ''),Chromosome_m,BAM_PATH)
                    TotalCov_m = eval(str(ELEPHANT.replace('-', ''))+"_CovDICT")[Chromosome_m]
                    
            try:
                SVcov = subprocess.check_output(ProcessToRun_SVcov, shell=True)
                CovProp = float(SVcov)/float(TotalCov)
                if float(TotalCov) > 0:
                    if CovProp > 1.8 and CovProp < 2.8:
                        SV_CoverageData.append("Homo")
                    elif CovProp < 1.8 and CovProp > 0.8:
                        SV_CoverageData.append("Het")
                    elif CovProp < 0.8 and CovProp > 0.3:
                        SV_CoverageData.append("Het*")
                    elif CovProp < 0.3:
                        SV_CoverageData.append("Deletion")    
                    elif CovProp > 2.8:
                        SV_CoverageData.append("Homo*")
                    else:
                        SV_CoverageData.append("OoR")
                else:
                    SV_CoverageData.append("NA")
            except subprocess.CalledProcessError as e:
                print(e.output)
                SV_CoverageData.append("Error")

            if SV_TYPE == "Rearrangements":
                try:
                    SVcov_m = subprocess.check_output(ProcessToRun_SVcov_m, shell=True)
                    CovProp_m = float(SVcov_m)/float(TotalCov_m)
                    if CovProp_m < CovProp_m:
                        PTally += 1
                    else:
                        MTally += 1
                except subprocess.CalledProcessError:
                    PTally += 1
        elif data_[5] != "Fail":
            SV_CoverageData.append("CovFail")       
        else:
            SV_CoverageData.append("Polarize")
        
        
    if MTally > PTally:
        RESwap = "True"
    else:
        RESwap = "False"

    if read_type == "Rearrangements":
        return SV_CoverageData, RESwap
    else:
        return SV_CoverageData
                         
        
        
def ElephantNamesToList(names_):
    names_ = names_.replace('[', '')
    names_ = names_.replace(']', '')
    names_ = names_.replace('\'', '')
    names_ = names_.replace(' ', '')
    names_ = names_.split(",")
    return names_



FINAL_DATA = []

AllElephantList = snakemake.params[2]

with open(SVData_Genotyped_Output, 'w') as f:
    write = csv.writer(f)
    
    with open(SVData_input, 'r') as file:
        line_reader = csv.reader(file, delimiter=',')
        
        for line in line_reader:
            if line[0] != "Chromosome":
                if SV_TYPE == "Rearrangements":
                    x,y = GetHeterozygosity(line, SV_TYPE)
                    write.writerow(line + [x]+[y]) #
                else:
                    x = GetHeterozygosity(line, SV_TYPE)
                    write.writerow(line + [x]) #
            else:
                if SV_TYPE == "Rearrangements":
                    write.writerow(line + ["Genotype", "SwapPnM"]) #
                else:
                    write.writerow(line + ["Genotype"])