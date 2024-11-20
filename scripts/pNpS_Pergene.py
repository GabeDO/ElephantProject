
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 14:34:53 2024

@author: goreill1
"""

import allel
import scipy
from collections import defaultdict
import gzip
import csv

#SV data import
#get list of geneIDs hit: [GeneID:1,NotOnGene,GeneID:3,Adj,GeneID:10]

path_vcf = str(snakemake.input[0]) #path to the vcf file
GTF = str(snakemake.input[1])
DNDS_path = str(snakemake.input[2])
TD_data = str(snakemake.input[3])
RE_data = str(snakemake.input[4])

OutputFile = str(snakemake.output[0])
#GTF data import
#make Dict from list of geneIDs containign exon pos: geneChromDict[GeneID:1] = [[1,4],[7,10]]

ChromList = []
with gzip.open(GTF, "rt", newline = '') as file: ##with gzip.open (GTF, "rt", newline = ") as file:
    rows = csv.reader (file, delimiter='\t')
    for row in rows:
        if len(row) != 1 and row[8].split("\"")[5][0] != "(":
            if row[0][0:2] != "NW":
                if row[2] == "gene": #CHECK IF ITS ONE OF THE TYPES I CARE
                    if "Gene_"+row[0].replace(".", "_") in locals (): # IF THE CHROMDICTONARY EXISTS
                        eval("Gene_"+row[0].replace(".", "_"))[row[8].split("\"")[5]] += [int(row[3]),int(row[4])]
                    else: # IF THE CHROMDICTONARY DOSEN'T EXISTS
                        print(row[0].replace(".", "_"))
                        ChromList.append(row[0])
                        exec("Gene_"+row[0].replace(".", "_") + "= defaultdict(list)") #MAKE A DICTONARY NAMED AFTER THE CHROMOSOME 
                        eval("Gene_"+row[0].replace(".", "_"))[row[8].split("\"")[5]] += [int(row[3]),int(row[4])] #ADD TO THE DICTONARY IN THE FORMAT: CHROM= {GENEID: START, END}


#SnNs data import

outputlist = []
with open(DNDS_path, 'r') as read_obj: 
    rows = csv.reader (read_obj, delimiter=',')
    for row in rows:
        if row[0] != "Chrom": # skip the header
            if row[0][0:2] != "NW":
                #if row[2] == "Sn" or row[2] == "Sn*" or row[2] == "4fS":
                if "DNDS_"+row[0].replace(".", "_") in locals (): # IF THE CHROMDICTONARY EXISTS
                    eval("DNDS_"+row[0].replace(".", "_"))[int(row[1])] = row[2]
                else: # IF THE CHROMDICTONARY DOSEN'T EXISTS
                    # just prints out the Chrom name. If all is working, each chrom is printed once
                    exec("DNDS_"+row[0].replace(".", "_") + "= defaultdict(list)") #MAKE A DICTONARY NAMED AFTER THE CHROMOSOME 
                    eval("DNDS_"+row[0].replace(".", "_"))[int(row[1])] = row[2] #ADD TO THE DICTONARY IN THE FORMAT: Sn_CHROM= {POS: Sn}


for Chromosome in ChromList:# iterate over chromososmes found in the gtf
    print("Calculating pNpS on chromosome: ",Chromosome)
    
    data_vcf = allel.read_vcf(path_vcf, region=Chromosome) #sci-kits thing to read a vcf into data, set the chom its looking at
    
    for i in range(len(data_vcf['calldata/GT'])): #this iterates over each index/row of the vcf file 
        try:
            found = False #set up a variable that'll change to true if it finds anything
            for key, values in eval("Gene_"+Chromosome.replace(".", "_")).items(): #iterate through the list of genes in the specified chromosome
                if int(data_vcf['variants/POS'][i]) > int(values[0]) and int(data_vcf['variants/POS'][i]) < int(values[1]): #if the current position is on a gene, then:
                    if "GeneID" in data_vcf: #its not a defult dictionary, so I gotta check if they key exists yet
                        found = True #set foudn to ture if something is found
                        data_vcf["GeneID"] += [key] #add the gene ID "[key]" to the new VCF "GeneID" key
                        break #if it finds something, we gucci. break the loop
                    else:
                        found = True
                        data_vcf["GeneID"] = [key]
                        break
            if found == False: #if nothing is found, then that position isn't on a gene. 
                if "GeneID" in data_vcf:
                    data_vcf["GeneID"] += ["Not Gene"]
                else:
                    data_vcf["GeneID"] = ["Not Gene"]
    
    
        except NameError:
            pass 
        
        try:
            if int(data_vcf['variants/POS'][i]) in eval("DNDS_"+(data_vcf['variants/CHROM'][i]).replace(".", "_")):
                if "DNDS" in data_vcf:
                    data_vcf["DNDS"] += [eval("DNDS_"+(data_vcf['variants/CHROM'][i]).replace(".", "_"))[data_vcf['variants/POS'][i]]]
                else:
                    data_vcf["DNDS"] = [eval("DNDS_"+(data_vcf['variants/CHROM'][i]).replace(".", "_"))[data_vcf['variants/POS'][i]]]
            else:
                if "DNDS" in data_vcf:
                    data_vcf["DNDS"] += ["NA"]
                else:
                    data_vcf["DNDS"] = ["NA"]
        except NameError:
            pass
        
        
    GeneList = list(set(data_vcf["GeneID"]))
    for i in GeneList:
        poslistNS = []
        poslistSn = []
        
        for geneID in range(len(data_vcf["GeneID"])):
            try:
                if data_vcf["GeneID"][geneID] == i and (data_vcf['DNDS'][geneID] == "Sn" or data_vcf['DNDS'][geneID] == "Sn*" or data_vcf['DNDS'][geneID] == "4fS"):
                    poslistSn.append(data_vcf['variants/POS'][geneID])
                    
                if data_vcf["GeneID"][geneID] == i and data_vcf['DNDS'][geneID] == "NS":
                    poslistNS.append(data_vcf['variants/POS'][geneID])
            except KeyError:
                pass
        
        if len(poslistNS) > 1 and  len(poslistSn) > 1:
                
            GT = allel.GenotypeArray(data_vcf['calldata/GT'])
            ac = GT.count_alleles()
            pN = allel.sequence_diversity(poslistNS,ac)
            pS = allel.sequence_diversity(poslistSn,ac)
            outputlist.append([Chromosome,i,pN,pS,pN/pS,len(poslistNS),len(poslistSn)])
        elif len(poslistNS) > 1:
            GT = allel.GenotypeArray(data_vcf['calldata/GT'])
            ac = GT.count_alleles()
            pN = allel.sequence_diversity(poslistNS,ac)
            outputlist.append([Chromosome,i,pN,'NA','NA',len(poslistNS),len(poslistSn)])
        elif len(poslistSn) > 1:
            GT = allel.GenotypeArray(data_vcf['calldata/GT'])
            ac = GT.count_alleles()
            pS = allel.sequence_diversity(poslistSn,ac)
            outputlist.append([Chromosome,i,'NA',pS,'NA',len(poslistNS),len(poslistSn)])
        else:
            outputlist.append([Chromosome,i,'NA','NA','NA',len(poslistNS),len(poslistSn)])

GeneIDList_TD = defaultdict(int)
with open(TD_data, "rt", newline = '') as file:
    rows = csv.reader(file, delimiter=',')
    data = list(rows)
    
    for row in data:
        GeneIDList_TD[row[8]] += 1

GeneIDList_RE = defaultdict(int)
with open(RE_data, "rt", newline = '') as file:
    rows = csv.reader(file, delimiter=',')
    data = list(rows)
    
    for row in data:
        GeneIDList_RE[row[13]] += 1

for r in outputlist:
    if r[1] in GeneIDList_TD:
        r.append(1)
    else:
        r.append(0)
    if r[1] in GeneIDList_RE:
        r.append(1)
    else:
        r.append(0)


with open(OutputFile, "w") as f:
    writer = csv.writer(f)
    writer.writerow(["Chromosome","GeneID","pN","pS","pNpS","Ns_count","Sn_count","TD","RE"]) 
    writer.writerows(outputlist)   
