#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  3 19:24:16 2025

@author: Gabe
"""

import gzip
import sys
import re
import statistics
import os
import csv
import numpy as np
from collections import defaultdict

def open_smart(filename, mode = "rt"):
    if filename.endswith(".gz"):
        return(gzip.open(filename, mode, newline = ''))
    else:
        return(open(filename, mode, newline = ''))


def CheckSV(SVStart, SVEnd, InputStart, InputEnd, type_, Start_, End_):
    Start = Start_
    End = End_
    if SVStart > InputStart and SVStart < InputEnd:
        Start = "SV starts at an " + str(type_)
        #print("COCA COLA?!")
        
    if SVEnd > InputStart and SVEnd < InputEnd:
        End = "SV ends at an " + str(type_)
        #print("YIPPE!")
        
    return Start, End


GTF = str(snakemake.input[0])

with open_smart(GTF, "rt") as file: ##with gzip.open (GTF, "rt", newline = ") as file:
    rows = csv.reader (file, delimiter='\t')
    for row in rows:
        try:
            if len(row) != 1 and row[8].split("\"")[5][0] != "(":
                if row[2] == "gene" or row[2] == "exon" or row[2] == "start_codon" or row[2] == "stop_codon": #CHECK IF ITS ONE OF THE TYPES I CARE
        
                    ## I have to use 'row [0].replace(".", "_")' to turn the "." in
                    ##the chrom names to "_" because python dosen't like when variables have "." in their names
                    if row[0].replace(".", "_") in locals (): # IF THE CHROMDICTONARY EXISTS
                        #ADD TO THE DICTONARY IN THE FORMAT: CHROM= {GENEID: TYPE, START, END, LENGTH} 
                        eval(row[0].replace(".", "_"))[row[8].split("\"")[5]] += [ [row[2], int (row [3]), int(row[4]), int(row[4])-int(row[3])] ]
        
                    else: # IF THE CHROMDICTONARY DOSEN'T EXISTS
                        print(row[0].replace(".", "_")) # just prints out the Chrom name. If all is working, each chrom is printed once
                        exec(row[0].replace(".", "_") + "= defaultdict(list)") #MAKE A DICTONARY NAMED AFTER THE CHROMOSOME 
                        #ADD TO THE DICTONARY IN THE FORMAT: CHROM= {GENEID: TYPE, START, END, LENGTH}
                        eval(row[0].replace(".", "_"))[row[8].split("\"")[5]] += [ [row[2], int (row [3]), int(row[4]), int(row[4])-int(row[3])] ]
        except IndexError:
            if len(row) != 1 and row[8].split(";")[2][0] != "(":
                if row[2] == "gene" or row[2] == "exon" or row[2] == "start_codon" or row[2] == "stop_codon": #CHECK IF ITS ONE OF THE TYPES I CARE
        
                    ## I have to use 'row [0].replace(".", "_")' to turn the "." in
                    ##the chrom names to "_" because python dosen't like when variables have "." in their names
                    if row[0].replace(".", "_") in locals (): # IF THE CHROMDICTONARY EXISTS
                        #ADD TO THE DICTONARY IN THE FORMAT: CHROM= {GENEID: TYPE, START, END, LENGTH} 
                        eval(row[0].replace(".", "_"))[row[8].split(";")[2]] += [ [row[2], int (row [3]), int(row[4]), int(row[4])-int(row[3])] ]
        
                    else: # IF THE CHROMDICTONARY DOSEN'T EXISTS
                        print(row[0].replace(".", "_")) # just prints out the Chrom name. If all is working, each chrom is printed once
                        exec(row[0].replace(".", "_") + "= defaultdict(list)") #MAKE A DICTONARY NAMED AFTER THE CHROMOSOME 
                        #ADD TO THE DICTONARY IN THE FORMAT: CHROM= {GENEID: TYPE, START, END, LENGTH}
                        eval(row[0].replace(".", "_"))[row[8].split(";")[2]] += [ [row[2], int (row [3]), int(row[4]), int(row[4])-int(row[3])] ]


CSVfile = str(snakemake.input[1])
SV_TYPE = str(snakemake.params[0])

SVList = []
SVList_mate =[]


ResultsDic = defaultdict(list)
ResultsDic_m = defaultdict(list)
GeneHitCounter = defaultdict(int)


Ticker = 0
with open(CSVfile, 'r') as file:
    rows = csv.reader(file)
    for SV in rows:
        if SV[0] == "Chromosome":
            if SV_TYPE == 'Rearrangements':
                OutputHeader = [SV+['GeneID','GeneLength','ExonNumber','ExonsHit','StartCodonHit','StopCodonHit','MutationInfo']+['GeneID_m','GeneLength_m','ExonNumber_m','ExonsHit_m','StartCodonHit_m','StopCodonHi_m','MutationInfo_m']]
            else:
                OutputHeader = [SV+['GeneID','GeneLength','ExonNumber','ExonsHit','StartCodonHit','StopCodonHit','MutationInfo']]
        elif SV[0] == str(snakemake.params[1]): #check if the pair chrom is the one snakemake is looking at atm (to split jobs by chr in snakemake)
            
            #make an ID for each SV: "chom_start-end"
            Result,Result_m = ['Not on Gene',0,0,0,0,0,'None'], ['Not on Gene',0,0,0,0,0,'None'] 
            if SV_TYPE == 'Rearrangements':
                SVID = SV[0]+"-"+str(SV[5])+"-"+str(SV[6])
                #Track where the SV starts abd ends
                SVStart = int(SV[5])
                SVEnd = int(SV[6])
                
                
                SVStart_m = int(SV[7])
                SVEnd_m = int(SV[8])
                
            elif SV_TYPE == 'TandemDups':
                SVID = SV[0]+"-"+str(SV[3])+"-"+str(SV[4])
                #Track where the SV starts abd ends
                SVStart = int(SV[3])
                SVEnd = int(SV[4])

                
            
            
            #check if the SV chromosome is found in the GTF. If it isn't, skip it
            
            if SV[0].replace(".", "_") in locals ():
                #now iterate through the dictionary of genes in that chrom
                for key, values in eval(SV[0].replace(".", "_")).items():
                    #get the gene length and save it for later
                    GeneLength = values[0][-1]
                    ExonNumber = len(values)
                
                    #check to see if the SV actually overlaps on the gene at all. If not, skip this gene.
                    if len(range(max(SVStart,values[0][1]), min(SVEnd,values[0][2]))) > 0:
                        
                        GeneHitCounter[key] += 1
            
                        #if it does overlap, see if the SV covers the whole gene. If so, just mark it as "Whole Gene" and call it a day
                        if SVStart < values[0][1] and SVEnd > values[0][2] and values[0][0] == "gene":
                            Result = [key,GeneLength,ExonNumber,'NA','NA','NA',"Whole Gene"] 
                            ResultsDic[SVID] = SV + Result
                        #otherwise, iterate through all the exons/codons/etc
                        else:
                            #set up a little value table
                            Hits = [key,GeneLength,ExonNumber,0,0,0]
                            Start_Info = "SV starts at an Intron"
                            End_Info = "SV ends at an Intron"
                            for g in values:
                                #add to the value table based on the type of data
                                #Check if its overlapping on an exon, and if so add to the results table
                                if len(range(max(SVStart,g[1]), min(SVEnd,g[2])+1)) > 0 and g[0] in ["exon", "start_codon","stop_codon"]:
                                    if g[0] == "exon":
                                        Hits[3] += 1
                                        Start_Info, End_Info = CheckSV(SVStart, SVEnd, g[1], g[2], g[0],Start_Info,End_Info)
                                    elif g[0] == "start_codon":
                                        Hits[4] += 1
                                        Start_Info, End_Info = CheckSV(SVStart, SVEnd, g[1], g[2], g[0],Start_Info,End_Info)
                                    elif g[0] == "stop_codon":
                                        Hits[5] += 1
                                        Start_Info, End_Info = CheckSV(SVStart, SVEnd, g[1], g[2], g[0],Start_Info,End_Info)
                                        
                                    
                            if sum(Hits[3:]) > 0: #Check if it actually hiyt anything. No intron-intron SVs unless theres an exon inbetween
                                Result = Hits + [Start_Info + " " + End_Info]
                            else: #if not, it just hits an intron
                                Result = [key,GeneLength,ExonNumber,0,0,0,"Intronic"]
                                
                    #A quick elseif to see if its in the range of the gene if its extended by 10,000 - checking if its adjacent.
                    elif len(range(max(SVStart,values[0][1]-10000), min(SVEnd,values[0][2]+10000))) > 0:
                        Result = [key,GeneLength,ExonNumber,0,0,0,"Adjacent"]
                        
                        
            if SV_TYPE == 'Rearrangements' and SV[2] == "=":
                if SV[0].replace(".", "_") in locals ():
                    for key, values in eval(SV[0].replace(".", "_")).items():
                        #get the gene length and save it for later
                        GeneLength = values[0][-1]
                        ExonNumber = len(values)
                        if len(range(max(SVStart_m,values[0][1]), min(SVEnd_m,values[0][2]))) > 0:
                            #if it does overlap, see if the SV covers the whole gene. If so, just mark it as "Whole Gene" and call it a day
                            if SVStart_m < values[0][1] and SVEnd_m > values[0][2] and values[0][0] == "gene":
                                Result_m = [key,GeneLength,ExonNumber,'NA','NA','NA',"Whole Gene"]
                                ResultsDic_m[SVID] += Result_m
                            #otherwise, iterate through all the exons/codons/etc
                            else:
                                #set up a little value table
                                Hits = [key,GeneLength,ExonNumber,0,0,0]
                                Start_Info = "SV starts at an Intron"
                                End_Info = "SV ends at an Intron"
                                for g in values:
                                    #add to the value table based on the type of data
                                    #Check if its overlapping on an exon, and if so add to the results table
                                    if len(range(max(SVStart_m,g[1]), min(SVEnd_m,g[2])+1)) > 0 and g[0] in ["exon", "start_codon","stop_codon"]:
                                        if g[0] == "exon":
                                            Hits[3] += 1
                                            Start_Info, End_Info = CheckSV(SVStart_m, SVEnd_m, g[1], g[2], g[0],Start_Info,End_Info)
                                        elif g[0] == "start_codon":
                                            Hits[4] += 1
                                            Start_Info, End_Info = CheckSV(SVStart_m, SVEnd_m, g[1], g[2], g[0],Start_Info,End_Info)
                                        elif g[0] == "stop_codon":
                                            Hits[5] += 1
                                            Start_Info, End_Info = CheckSV(SVStart_m, SVEnd_m, g[1], g[2], g[0],Start_Info,End_Info)
                                            
                                        
                                if sum(Hits[3:]) > 0: #Check if it actually hiyt anything. No intron-intron SVs unless theres an exon inbetween
                                    Result_m = Hits + [Start_Info + " " + End_Info]
                                else: #if not, it just hits an intron
                                    Result_m = [key,GeneLength,ExonNumber,0,0,0,"Intronic"]
                                    
                        #A quick elseif to see if its in the range of the gene if its extended by 10,000 - checking if its adjacent.
                        elif len(range(max(SVStart_m,values[0][1]-10000), min(SVEnd_m,values[0][2]+10000))) > 0:
                            Result_m = [key,GeneLength,ExonNumber,0,0,0,"Adjacent"]
                    
            elif SV_TYPE == 'Rearrangements' and SV[2] != "=":
                if SV[2].replace(".", "_") in locals ():
                    for key, values in eval(SV[2].replace(".", "_")).items():
                        #get the gene length and save it for later
                        GeneLength = values[0][-1]
                        ExonNumber = len(values)
                        if len(range(max(SVStart_m,values[0][1]), min(SVEnd_m,values[0][2]))) > 0:
                            #if it does overlap, see if the SV covers the whole gene. If so, just mark it as "Whole Gene" and call it a day
                            if SVStart_m < values[0][1] and SVEnd_m > values[0][2] and values[0][0] == "gene":
                                Result_m = [key,GeneLength,ExonNumber,'NA','NA','NA',"Whole Gene"]
                                ResultsDic_m[SVID] += Result_m
                            #otherwise, iterate through all the exons/codons/etc
                            else:
                                #set up a little value table
                                Hits = [key,GeneLength,ExonNumber,0,0,0]
                                Start_Info = "SV starts at an Intron"
                                End_Info = "SV ends at an Intron"
                                for g in values:
                                    #add to the value table based on the type of data
                                    #Check if its overlapping on an exon, and if so add to the results table
                                    if len(range(max(SVStart_m,g[1]), min(SVEnd_m,g[2])+1)) > 0 and g[0] in ["exon", "start_codon","stop_codon"]:
                                        if g[0] == "exon":
                                            Hits[3] += 1
                                            Start_Info, End_Info = CheckSV(SVStart_m, SVEnd_m, g[1], g[2], g[0],Start_Info,End_Info)
                                        elif g[0] == "start_codon":
                                            Hits[4] += 1
                                            Start_Info, End_Info = CheckSV(SVStart_m, SVEnd_m, g[1], g[2], g[0],Start_Info,End_Info)
                                        elif g[0] == "stop_codon":
                                            Hits[5] += 1
                                            Start_Info, End_Info = CheckSV(SVStart_m, SVEnd_m, g[1], g[2], g[0],Start_Info,End_Info)
                                            
                                        
                                if sum(Hits[3:]) > 0: #Check if it actually hiyt anything. No intron-intron SVs unless theres an exon inbetween
                                    Result_m = Hits + [Start_Info + " " + End_Info]
                                else: #if not, it just hits an intron
                                    Result_m = [key,GeneLength,ExonNumber,0,0,0,"Intronic"]
                        
                    
                    
                    
            if SV_TYPE == 'Rearrangements':        
                ResultsDic[SVID] = [SV+Result+Result_m]
            else:
                ResultsDic[SVID] = [SV+Result]
                        
            if ResultsDic[SVID] == []:
                if SV_TYPE == 'Rearrangements':
                    ResultsDic[SVID] = [SV+Result+Result_m]
                else:
                    ResultsDic[SVID] = [SV+Result]
    

    
Outputfile = str(snakemake.output)

with open(Outputfile, 'w', newline='') as f: 
    writer = csv.writer(f, delimiter=',')
    writer.writerows(OutputHeader)
    
    for key, values in ResultsDic.items():
        writer.writerows(ResultsDic[key])

        
