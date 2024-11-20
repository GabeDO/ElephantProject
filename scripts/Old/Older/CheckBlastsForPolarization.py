import gzip
import sys
import re
import statistics
import os
import csv
from collections import defaultdict


pairInput = str(snakemake.input[0])
pairBlast = defaultdict(list)
with open(pairInput, "rt") as csv_file:
    csv_reader = csv.reader(csv_file, delimiter='\t')
    for row in csv_reader:
        if row[0] in pairBlast:
            pairBlast[row[0]] += [[row[1],row[8],row[9]]]
        else:
            pairBlast[row[0]] = [[row[1],row[8],row[9]]]
#Clear some memory:
del(csv_reader)


mateInput = str(snakemake.input[1])
mateBlast = defaultdict(list)
with open(mateInput, "rt") as csv_file:
    csv_reader = csv.reader(csv_file, delimiter='\t')
    for row in csv_reader:
        if row[0] in mateBlast:
            mateBlast[row[0]] += [[row[1],row[8],row[9]]]
        else:
            mateBlast[row[0]] = [[row[1],row[8],row[9]]]
#Clear some memory:
del(csv_reader)

PolarizedList = defaultdict(list)

DistanceThreshold = 1000
OverlapThreshold = 1000
#for each mutation blast found

def CheckValidBlast(minA, maxA, minB, maxB, DisThresh, OvrThresh):
    checkvalue = 0
    if len(range(max(minA, minB), (min(maxA, maxB)+DisThresh))) > 0:
        checkvalue += 1
    if len(range(max(minA, minB), (min(maxA, maxB)-OvrThresh))) > 0:
        checkvalue += 1
    if not (minA > minB and maxA < maxB):
        checkvalue += 1
    if not (minB > minA and maxB < maxA):
        checkvalue += 1

    if checkvalue == 4:
        return True
    else:
        return False


for i in pairBlast:
    #check if it has its mate in there
    if i in mateBlast and i not in PolarizedList:
        #then loop through...
        for p_hits in pairBlast[i]:
            if i in PolarizedList:
                break
            
            
            ##Check to see if the 'MateTooFarAway' condition has been met (see below) - only works if p_hits and m_hits are sorted by [1]
            #if MateTooFarAway == True:
            #    MateTooFarAway = False
            #    break
            
            
            #each blast hit in them. 
            for m_hits in mateBlast[i]:
                if i in PolarizedList:
                    break
                
                #checking to see if the search has passed - only works if p_hits and m_hits are sorted by [1]
                #if (m_hits[2] - p_hits[1] > DistanceThreshold):
                #    MateTooFarAway = True
                #    break
                
                #if the hits are on the same chrom
                if p_hits[0] == m_hits[0]:
                    #check if they're close to each other (within 100 bp)
                    if CheckValidBlast(int(p_hits[1]), int(p_hits[2]),int(m_hits[1]),int(m_hits[2]), DistanceThreshold, OverlapThreshold) == True:
                        PolarizedList[i] = "GG Got em"
#remove bigs lists once i'm done witht them to save memory:
del(pairBlast)
del(mateBlast)





NonPolarizedData = str(snakemake.input[2])
Newdata = []
with open(NonPolarizedData, "rt") as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
        Mutation = row[0]+":"+str(row[5])+"-"+str(row[6])+"_"+row[2]+":"+str(row[7])+"-"+str(row[8])
        Mutation = Mutation.replace(" ", "")
        if row[0] == "Chromosome":
            row.insert(9,"Polarization")
            Newdata.append(row)
        elif Mutation in PolarizedList:
            row.insert(9,"Polar")
            Newdata.append(row)
        else:
            row.insert(9,"NotPolar")
            Newdata.append(row)
#Clear some memory:
del(csv_reader)

with open(r'Mutations_to_polarize.csv', 'w') as fp:
    for item in Newdata:
        # write each item on a new line
        if item[0] == "Chromosome":
            fp.write("Chromosome,ReadBinPosition,MateChromosome,MateBinPosition,NumberOfElephants,PairMin,PairMax,MateMin,MateMax,Polarized,CoverageList,ElephantNames\n")
        else:
            fp.write("{},{},{},{},{},{},{},{},{},{},{},{}\n".format(item[0],int(item[1]),item[2],int(item[3]),int(item[4]),int(item[5]),int(item[6]),int(item[7]),int(item[8]),item[9],item[10],item[11]))

#[0]+":"+str(5)+"-"+str(6)+"_"[2]+":"+str([7])+"-"str([8])
#['NC_064819.1', '101770200', 'NC_064839.1', '36018500', '2', ' 101770696', ' 101770849', ' 36019312', ' 36019463', 'NotPolar', "['Female-127_S10', 'Female-15_S1']"]
#NC_064819.1:101770696-101770849_NC_064839.1:36019312-36019463
#NC_064819.1,101770200,NC_064839.1,36018500,2, 101770696, 101770849, 36019312, 36019463,"['Female-127_S10', 'Female-15_S1']"