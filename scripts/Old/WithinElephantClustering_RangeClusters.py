import gzip
import sys
import re
import statistics
import os
from collections import defaultdict



Dict = defaultdict(list)
binsize = 700
ReadLen = 151
OutputList = []

PreviousReadPosition = 0
PreviousMatePosition = 0
INPUT = open(str(snakemake.input), 'r')

PreviousBin = 0
PreviousMateBin = 0
ReadsInMutation = []


def checkBin(Value, prevValue):
    Data = Value.split(',')
    PrevData = prevValue.split(',')
    ErrorList = []
    
    if Value[0] == prevValue[0]:
        if Data[1] < PrevData[2]:
            ErrorList.append([Data,"Reads"])
        if Data[4] < PrevData[5]:
            ErrorList.append([Data, "Mates"])
        with open("BinErrorsRead.txt", 'a') as Err:  
            if len(ErrorList) > 1:
                for i in ErrorList:
                    Err.write(i)
            else:
                pass
    return


for line in INPUT:
    
    A=line.split(' ')
    Chromosome = A[0]
    ReadPosition = int(A[1])
    MateChromosome = A[2]
    MatePosition = int(A[3])
    Elephant = A[-1]
    if len(ReadsInMutation) == 0:
        ReadsInMutation = [ReadPosition]
        MateReadsInMutation = [MatePosition]

    elif ReadPosition <= PreviousReadPosition + (350 + ReadLen) and PreviousBin != (min(ReadsInMutation)-350) - (min(ReadsInMutation)-350) % (binsize):
        ReadsInMutation.append(ReadPosition)
        MateReadsInMutation.append(MatePosition)
    else:       
        MinReadbinPosition = (min(ReadsInMutation)-350) - (min(ReadsInMutation)-350) % (binsize)
        MinMateBinPosition = (min(MateReadsInMutation)) - (min(MateReadsInMutation)-350) % (binsize)
        
        MaxReadbinPosition = (max(ReadsInMutation)+ReadLen-350) - (max(ReadsInMutation)+ReadLen-350) % (binsize)
        MaxMateBinPosition = (max(MateReadsInMutation)+ReadLen-350) - (max(MateReadsInMutation)+ReadLen-350) % (binsize)

        #PUT ANOTHER FUCNTION HERE TO SORT THEN ITERATE THROUGH MateReadsInMutation THEN SPLIT THE BINS UP AT AREAS THEY AREN'T WITHIN 350 OF EACHOTHER

        BinKey = f"{Chromosome},{MinReadbinPosition},{MaxReadbinPosition},{MateChromosome},{MinMateBinPosition},{MaxMateBinPosition}"    
        ReadsInMutation = [ReadPosition]
        MateReadsInMutation = [MatePosition]

        Dict[BinKey] += ReadsInMutation

    PreviousBin = (min(ReadsInMutation)-350) - (min(ReadsInMutation)-350) % (binsize)
    PreviousReadPosition = int(A[1])
    PreviousMatePosition = int(A[3])

with open(str(snakemake.output), 'w') as f:
    prevvalue = "ChrYEET,0,0,ChrMEET,0,0"

    for key, value in Dict.items():
        if len(value) >= 3:
            f.write('%s,%s,%s,%s' % (key, value, Elephant))
            checkBin(value, prevvalue)        
            prevvalue = value
        else:
            pass
