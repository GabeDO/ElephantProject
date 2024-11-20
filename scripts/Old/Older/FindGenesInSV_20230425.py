import gzip
import sys
import re
import statistics
import os
import csv
import numpy as np
from collections import defaultdict

#How close a rearrangement needs to be to a gene for me to consider it importants to expression
##OverlapThreshhold = int(snakemake.params[0]) + 1

##GeneDic = defaultdict(list)

#GTF = '/Users/goreill1/Desktop/GeneIdentification/GCF_024166365.1_mEleMax1_primary_haplotype_genomic.gtf.gz'
##GTF = str(snakemake.input[0])
##with gzip.open(GTF,"rt", newline = '') as file:
##    rows = csv.reader(file, delimiter='\t')
##    for row in rows:
##        if len(row) != 1:
##            if row[2] == "gene":
                
##                GeneDic[row[8].split("\"")[5]] = [row[0],int(row[3]),int(row[4])]
                

#CSVfile = '/Users/goreill1/Desktop/GeneIdentification/Rearrangements_all_100kDistance.csv'
##CSVfile = str(snakemake.input[1])
##with open(CSVfile, 'r') as csv_file:
##    rows = csv.reader(csv_file)
##    for row in rows:
##        143898
##        if row[0] == 'Chromosome':
##            pass
##        else:
##            ReadRange = [int(row[5]),int(row[6])]
##            PairRange = [int(row[7]),int(row[8])]
            
##            Chromosome = row[0]
##            for key, values in GeneDic.items():
##                #CheckOverlap < 1 means overlap
##                if len(range(max(ReadRange[0], values[1]), min(ReadRange[-1], values[2])+OverlapThreshhold)) > 0:
##                    #if len three, it hasnt had its rearrangment counter on it yet, so I'll odd one
##                    if len(values) == 3:
##                        GeneDic[key] += [[1,0]]
##                    #otherise, it already does have one, so I'll just add to the pair counter
##                    else:
##                        GeneDic[key][-1][0] += 1
##                if len(range(max(PairRange[0], values[1]), min(PairRange[-1], values[2])+OverlapThreshhold)) > 0:
##                    #same thing for the mate reads
##                    if len(values) == 3:
##                        GeneDic[key] += [[0,1]]
##                    else:
##                        GeneDic[key][-1][1] += 1
                

##Outputfile = str(snakemake.output)
##with open(Outputfile, 'w') as f: 
##    for key, value in GeneDic.items():
##        if len(value) != 3:
##            f.write('%s,%s,%s\n' % (key, value[-1], sum(value[-1])))



import gzip
import sys
import re
import statistics
import os
import csv
import numpy as np
from collections import defaultdict

Gen = 0
Rea = 0

GeneDic = defaultdict(list)

GTF = str(snakemake.input[0])
GeneList = []
with gzip.open(GTF,"rt", newline = '') as file:
    rows = csv.reader(file, delimiter='\t')
    for row in rows:
        if len(row) != 1:
            if row[2] == "gene" or row[2] == "exon" or row[2] == "start_codon" or row[2] == "stop_codon":
                GeneList.append([row[2],row[8].split("\"")[5],int(row[3]), int(row[4]), row[0], int(row[4])-int(row[3])])

#Key=GeneID: Start, End, [[Exon1 Start, Exon1 End],[Exon2 Start, Exon2 End]...]

print("exons found:")
print(len(GeneList))


CSVfile = str(snakemake.input[1])
SVList = []
with open(CSVfile, 'r') as file:
    rows = csv.reader(file)
    for row in rows:
        if row[0] != "Chromosome":
            SVList.append([[int(row[3]),int(row[4])], row[0]])
            #SV list items: StartPos, EndPos, Chrom

print("\nSVs found:")
print(len(SVList))

def SearchInputFile(Input, geneMin, geneMax, geneChrom, OverlapThreshhold, TYPE):
    GeneHits = 0
    WholeGene = 0
    PartialGene = 0
    DupList = []
    for row_ in Input:
        if row_[0] == 'Chromosome':
            pass
        else:
            ReadRange = row_[0]
            Chromosome = row_[1]

            #CheckOverlap < 1 means overlap
            if len(range(max(ReadRange[0], geneMin), min(ReadRange[-1], geneMax)+OverlapThreshhold)) > 0 and Chromosome == geneChrom:
                GeneHits += 1
                DupList.append(ReadRange)
                if TYPE == "gene":
                    if ReadRange[0] < geneMin and ReadRange[-1] > geneMax and Chromosome == geneChrom:
                        WholeGene += 1
                    else:
                        PartialGene += 1
                #Look for TD breakpoints in Exons
                elif TYPE == "exon":
                    StartBreakpoints = 0
                    EndBreakpoints = 1
                    if ReadRange[0] > geneMin and ReadRange[0] < geneMax:
                        StartBreakpoints += 1
                    if ReadRange[0] > geneMin and ReadRange[0] < geneMax:
                        EndBreakpoints += 1


    return [int(GeneHits), int(WholeGene), int(PartialGene), DupList]


for gene in GeneList:  
    GeneID = gene[1]
    if gene[0] == "gene":
        value = SearchInputFile(SVList, gene[2], gene[3], gene[4], 1, True)
        if GeneID in GeneDic and value[0] != 0:
            GeneDic[GeneID][0] += value[0]
            GeneDic[GeneID][4] += value[1]
            GeneDic[GeneID][5] += value[2]
            GeneDic[GeneID][5] += value[3]            
        elif value[0] != 0:
            GeneDic[GeneID] = [value[0],0,0,0,value[1],value[2], gene[-1],0,value[3]]
    elif gene[0] == "exon":
        if GeneID in GeneDic:
            GeneDic[GeneID][7] += 1
        else:
            GeneDic[GeneID] = [0,0,0,0,0,0, gene[-1],1,[]]
        
        value = SearchInputFile(SVList, gene[2], gene[3], gene[4], 1, False)
        if GeneID in GeneDic and value[0] != 0:
            GeneDic[GeneID][1] += value[0]
        elif value[0] != 0:
            GeneDic[GeneID] = [0,value[0],0,0,0,0, gene[-1],0,[]]
    elif gene[0] == "start_codon":
        value = SearchInputFile(SVList, gene[2], gene[3], gene[4], 1, False)
        if GeneID in GeneDic and value[0] != 0:
            GeneDic[GeneID][2] += value[0]
        elif value[0] != 0:
            GeneDic[GeneID] = [0,0,value[0],0,0,0, gene[-1],0,[]]
    elif gene[0] == "stop_codon":
        value = SearchInputFile(SVList, gene[2], gene[3], gene[4], 1, False)
        if GeneID in GeneDic and value[0] != 0:
            GeneDic[GeneID][3] += value[0]
        elif value[0] != 0:
            GeneDic[GeneID] = [0,0,0,value[0],0,0, gene[-1],0,[]]

Outputfile = str(snakemake.output)
with open(Outputfile, 'w') as f: 
    f.write('Gene,GeneHits,ExonHits,StartCodon,StopCodon,WholeHits,PartialHit,GeneSize,ExonNumber,TDList\n')
    for key, values in GeneDic.items():
        if len(values) != 3:
            f.write('%s,%s\n' % (key, values))