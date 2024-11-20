#'''
#Created on Wed May 24 11:54:47 2023
#@author: goreilli
#'''

import gzip
import sys
import re
import statistics
import os
import csv
import numpy as np
from collections import defaultdict

GTF = str(snakemake.input[0])
with gzip.open(GTF, "rt", newline = '') as file: ##with gzip.open (GTF, "rt", newline = ") as file:
    rows = csv.reader (file, delimiter='\t')
    for row in rows:
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



CSVfile = str(snakemake.input[1])
SV_TYPE = str(snakemake.params[0])

#with open (CSV file, 'r') as file:
SVList = []
SVList_mate =[]
with open(CSVfile, 'r') as file:
    rows = csv.reader(file)
    for row in rows:
        if row[0] != "Chromosome":
            if SV_TYPE == 'TandemDups':
                SVList.append([[int(row [3]), int(row [4])], row[0],row[-2],row[-3]])##This just makes a list of all the SVsI have.
                #the format is: [start, end] ,Chrom
            elif SV_TYPE == 'Rearrangements':
                if row[-1] == "True":
                    SVList_mate.append([[int(row [5]), int(row [6])], row[0],row[-2],row[-3]])
                    if row[2] == "=":
                        SVList.append([[int(row [7]), int(row [8])], row[0],row[-2],row[-3]])
                    else:
                        SVList.append([[int(row [7]), int(row [8])], row[2],row[-2],row[-3]])

                if row[-1] == "False":
                    SVList.append([[int(row [5]), int(row [6])], row[0],row[-2],row[-3]])
                    if row[2] == "=":
                        SVList_mate.append([[int(row [7]), int(row [8])], row[0],row[-2],row[-3]])
                    else:
                        SVList_mate.append([[int(row [7]), int(row [8])], row[2],row[-2],row[-3]])
                    

                 
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


#make dictonary for holding results
ResultsDic = defaultdict(list)
ResultsDic_m = defaultdict(list)
GeneHitCounter = defaultdict(int)

#iterate through my list of SV
for i in range(len(SVList)):
    
    #make an ID for each SV: "chom_start-end"
    
    SVID = SVList[i][1]+"-"+str(SVList[i][0][0])+"-"+str(SVList[i][0][1])
    #Track where the SV starts abd ends
    SVStart = SVList[i][0][0]
    SVEnd = SVList[i][0][1]
    SVGenoInfo = [SVList[-2],SVList[-1]]
    
    if SV_TYPE == 'Rearrangements':
        SVID_m = [SVList_mate[i][1]+"-"+str(SVList_mate[i][0][0])+"-"+str(SVList_mate[i][0][1])]
        SVStart_m = SVList_mate[i][0][0]
        SVEnd_m = SVList_mate[i][0][1]

    
    #check if the SV chromosome is found in the GTF. If it isn't, skip it
    if SVList[i][1].replace(".", "_") in locals ():
        #now iterate through the dictionary of genes in that chrom
        for key, values in eval(SVList[i][1].replace(".", "_")).items():
            #get the gene length and save it for later
            GeneLength = values[0][-1]
            ExonNumber = len(values)
        
            #check to see if the SV actually overlaps on the gene at all. If not, skip this gene.
            if len(range(max(SVStart,values[0][1]), min(SVEnd,values[0][2]))) > 0:
                
                GeneHitCounter[key] += 1

                #if it does overlap, see if the SV covers the whole gene. If so, just mark it as "Whole Gene" and call it a day
                if SVList[i][0][0] < values[0][1] and SVList[i][0][1] > values[0][2] and values[0][0] == "gene":
                    Result = [[key,GeneLength,ExonNumber,"Whole Gene"]] 
                    ResultsDic[i] += [Result,SVGenoInfo]
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
                        Result = [Hits] + [Start_Info, End_Info]
                        ResultsDic[i] += [Result,SVGenoInfo]
           #A quick elseif to see if its in the range of the gene if its extended by 10,000 - checking if its adjacent.
            elif len(range(max(SVStart,values[0][1]-10000), min(SVEnd,values[0][2]+10000))) > 0:
                Result = [[key,GeneLength,ExonNumber,"Adjacent"]]
                ResultsDic[i] += [Result,SVGenoInfo]
                
            if SV_TYPE == 'Rearrangements':
                if len(range(max(SVStart_m,values[0][1]), min(SVEnd_m,values[0][2]))) > 0:
                    #if it does overlap, see if the SV covers the whole gene. If so, just mark it as "Whole Gene" and call it a day
                    if SVStart_m < values[0][1] and SVEnd_m > values[0][2] and values[0][0] == "gene":
                        Result_m = [[key,GeneLength,ExonNumber,"Whole Gene"]] 
                        ResultsDic_m[i] += [Result_m,SVGenoInfo]
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
                            Result_m = [Hits] + [Start_Info, End_Info]
                            ResultsDic_m[i] += [Result_m,SVGenoInfo]
                #A quick elseif to see if its in the range of the gene if its extended by 10,000 - checking if its adjacent.
                elif len(range(max(SVStart_m,values[0][1]-10000), min(SVEnd_m,values[0][2]+10000))) > 0:
                    Result_m = [[key,GeneLength,ExonNumber,"Adjacent"]]
                    ResultsDic_m[i] += [Result_m,SVGenoInfo]
                
    if ResultsDic[i] == []:
        Result = ["Not on gene"]
        ResultsDic[i] = [Result,SVGenoInfo]

        

    
Outputfile = str(snakemake.output)
with open(Outputfile, 'w') as f: 
    
    f.write('ID,Chrom, GenesHit\n')
    
    for key, values in ResultsDic.items():

        f.write('%s,%s,%s\n' % (str(key)+"_p", SVList[key][1]+","+str(SVList[key][0][0])+","+str(SVList[key][0][1]) ,values))
        if SV_TYPE == 'Rearrangements':
            if ResultsDic_m[key] == []:
                ResultsDic_m[key] = ["Not on gene"]
            if SVList_mate[key][1] == "=":
                f.write('%s,%s,%s\n' % (str(key)+"_m", SVList[key][1]+","+str(SVList_mate[key][0][0])+","+str(SVList_mate[key][0][1]) ,ResultsDic_m[key]))
            else:
                f.write('%s,%s,%s\n' % (str(key)+"_m", SVList_mate[key][1]+","+str(SVList_mate[key][0][0])+","+str(SVList_mate[key][0][1]) ,ResultsDic_m[key]))

print('________________________________')
print(SV_TYPE)
print('________________________________')
for key, values in GeneHitCounter.items():
    print(str(key)+ ": "+ str(values))
print('________________________________')
print('________________________________')