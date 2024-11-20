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

                                                                          

VCFfile = str(snakemake.input[1])

#with open (CSV file, 'r') as file:
SNPList = []

with gzip.open(VCFfile, 'rt') as file:
    rows = csv.reader(file, delimiter='\t')
    for row in rows:
        if row[0][0] != "#":
                SNP_freq = row[7].split("AC=")[1].split(";")[0]
                SNPList.append([int(row[1]), row[0], SNP_freq])##This just makes a list of all the SVsI have. Chrom, Pos, Freq

         ########           



#make dictonary for holding results
ResultsDic = defaultdict(list)
GeneHitCounter = defaultdict(int)

#iterate through my list of SV
for i in range(len(SNPList)):
    
    #make an ID for each SV: "chom_start-end"
    
    SVID = SNPList[i][1]+"-"+str(SNPList[i][0])
    #Track where the SV starts abd ends
    SNPPos = SNPList[i][0]

    
    #check if the SNP chromosome is found in the GTF. If it isn't, skip it
    if SNPList[i][1].replace(".", "_") in locals ():
        #now iterate through the dictionary of genes in that chrom
        for key, values in eval(SNPList[i][1].replace(".", "_")).items():
            #get the gene length and save it for later
            GeneLength = values[0][-1]
            ExonNumber = len(values)
        
            #check to see if the SNP actually overlaps on the gene at all. If not, skip this gene.
            if SNPPos >= values[0][1] and SNPPos <= values[0][2]:                
                GeneHitCounter[key] += 1

                #set up a little value table
                Hits = [key,GeneLength,ExonNumber,0,0,0]
                Info = "SNP on an Intron"
                for g in values:
                    #add to the value table based on the type of data
                    #Check if its overlapping on an exon, and if so add to the results table
                    if SNPPos >= g[1] and SNPPos <= g[2] and g[0] in ["exon", "start_codon","stop_codon"]:
                        if g[0] == "exon":
                            Hits[3] += 1
                            Info = "SNP on an Exon"
                        elif g[0] == "start_codon":
                            Hits[4] += 1
                            Info  = "SNP on an Start Codon"
                        elif g[0] == "stop_codon":
                            Hits[5] += 1
                            Info = "SNP on an Stop Codon"
                            
                        
                if sum(Hits[3:]) > 0: #Check if it actually hiyt anything. No intron-intron SVs unless theres an exon inbetween
                    Result = [Hits] + [Info]
                    ResultsDic[i] += [Result]
                else:
                    Result = [Hits] + ["SNP on an Intron"]
                    ResultsDic[i] += [Result]

           #A quick elseif to see if its in the range of the gene if its extended by 10,000 - checking if its adjacent.
            elif SNPPos >= values[0][1]-10000 and SNPPos <= values[0][2]+10000 > 0:
                Result = [[key,GeneLength,ExonNumber,"Adjacent"]]
                ResultsDic[i] += [Result]
                
    if ResultsDic[i] == []:
        Result = ["Not on gene"]
        ResultsDic[i] = [Result]

        

    
Outputfile = str(snakemake.output)
with open(Outputfile, 'w') as f: 
    f.write('ID,Chrom,POS,Freq, GenesHit\n')
    for key, values in ResultsDic.items():
        f.write('%s,%s,%s,%s\n' % (str(key), SNPList[key][1]+","+str(SNPList[key][0]),SNPList[key][2],values,))