
'''
Created on Wed May 24 11:54:47 2023
@author: goreilli
'''

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
#with open (CSV file, 'r') as file:
SVList = []
with open(CSVfile, 'r') as file:
    rows = csv.reader(file)
    for row in rows:
        if row[0] != "Chromosome":
            SVList.append([[int(row [3]), int(row [4])], row[0]])##This just makes a list of all the SVsI have.
            #the format is: Chrom, start, end
            
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


#iterate through my list of SV
for i in SVList:
    
    #make an ID for each SV: "chom_start-end"
    SVID = i[1]+"_"+str(i[0][0])+"-"+str(i[0][1])
    #Track where the SV starts abd ends
    SVStart = i[0][0]
    SVEnd = i[0][1]
    
    #check if the SV chromosome is found in the GTF. If it isn't, skip it
    if i[1].replace(".", "_") in locals ():
    
        #now iterate through the dictionary of genes in that chrom
        for key, values in eval(i[1].replace(".", "_")).items():
            
            #get the gene length and save it for later
            GeneLength = values[0][-1]
            ExonNumber = len(values)
        
            #check to see if the SV actually overlaps on the gene at all. If not, skip this gene.
            if len(range(max(i[0][0],values[0][1]), min(i[0][1],values[0][2]))) > 0:
            
                #if it does overlap, see if the SV covers the whole gene. If so, just mark it as "Whole Gene" and call it a day
                if i[0][0] < values[0][1] and i[0][1] > values[0][2] and values[0][0] == "gene":
                    ResultsDic[SVID] += [[key,GeneLength,ExonNumber,"Whole Gene"]]
                
                #otherwise, iterate through all the exons/codons/etc
                else:
                    #set up a little value table
                    Hits = [key,GeneLength,ExonNumber,0,0,0]
                    Start_Info = "SV starts at an Intron"
                    End_Info = "SV ends at an Intron"
                    for g in values:
                        #add to the value table based on the type of data
                        #Check if its overlapping on an exon, and if so add to the results table
                        if len(range(max(i[0][0],g[1]), min(i[0][1],g[2])+1)) > 0 and g[0] in ["exon", "start_codon","stop_codon"]:
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
                        ResultsDic[SVID] += [Hits] + [Start_Info, End_Info] 
                            
            else:
                pass
    else:
        pass

Outputfile = str(snakemake.output)
with open(Outputfile, 'w') as f: 
    f.write('SV, GenesHit\n')
    for key, values in ResultsDic.items():
        f.write('%s,%s\n' % (key, values))