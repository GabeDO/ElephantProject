# REMOVE READLENGHTN FROM BINNING


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


def DataSummary(LIST):
    ListPair = []
    ListMate = []
    minPair = LIST[0][1]
    maxPair = LIST[0][1]+ReadLen
    minMate = LIST[0][3]
    maxMate = LIST[0][3]+ReadLen
    for i in LIST:
        ListPair.append(i[1])
        ListMate.append(i[3])
        if i[1] < minPair:
            minPair = i[1]
        elif i[1]+ReadLen > maxPair:
            maxPair = i[1]+ReadLen
        else:
            pass

        if i[3] < minMate:
            minMate = i[3]
        elif i[3]+ReadLen > maxMate:
            maxMate = i[3]+ReadLen
        else:
            pass

    return [[LIST[0][0],minPair,maxPair,(maxPair-minPair),LIST[0][2],minMate,maxMate,(maxMate-minMate)],[ListPair, ListMate]]

def checkBin(DataSum, key, PrevDataSum, prevkey, Elephant):
    with open("BinErrorsTEST.txt", 'a') as Err:  
        ReadError = ""
        MateError = ""
        if DataSum[0][0] == PrevDataSum[0][0] and DataSum[0][4] == PrevDataSum[0][4]:
            if PrevDataSum[0][1] <= DataSum[0][1] <= PrevDataSum[0][2] and PrevDataSum[0][5] <= DataSum[0][5] <= PrevDataSum[0][6]:
                if PrevDataSum[0][1] <= DataSum[0][1] <= DataSum[0][2] <= PrevDataSum[0][2]:
                    ReadError = "Read within previous"
                if PrevDataSum[0][5] <= DataSum[0][5] <= DataSum[0][6] <= PrevDataSum[0][6]:
                    MateError = "Mate within previous"
                if len(ReadError) + len(MateError) == 0:
                    Err.write(f"{key}:{DataSum[0]},Reads in bin:,{DataSum[1]} ||||| {prevkey}:{PrevDataSum[0]},prev reads in bin,{PrevDataSum[1]},{Elephant}")
                else:
                    pass
            else:
                pass

for line in INPUT:
    
    A=line.split(' ')
    Chromosome = A[0]
    ReadPosition = int(A[1])
    MateChromosome = A[2]
    MatePosition = int(A[3])
    Elephant = A[-1]
    
    
    #is read position is within 350+readlength of the last read, then:
    if ReadPosition <= PreviousReadPosition + (350 + ReadLen) and PreviousBin != (ReadPosition-350) - (ReadPosition-350) % (binsize + ReadLen):
        #put set its bin to be the same as the previous reads
        ReadbinPosition = PreviousBin
    else:
        #otherwise, keep going as normal and reset 'previous bin' values
        ReadbinPosition = (ReadPosition-350) - (ReadPosition-350) % (binsize + ReadLen)
        PreviousBin = (ReadPosition-350) - (ReadPosition-350) % (binsize + ReadLen)
        
    #Same thing, but for the mate reads 
    if MatePosition <= PreviousMatePosition + (350 + ReadLen) and PreviousMateBin != (MatePosition-350) - (MatePosition-350) % (binsize + ReadLen):
        MatebinPosition = PreviousMateBin
    else:
        PreviousMateBin = (MatePosition-350) - (MatePosition-350) % (binsize + ReadLen)
        MatebinPosition = (MatePosition-350) - (MatePosition-350) % (binsize + ReadLen)
    
    BinKey=f"{Chromosome},{ReadbinPosition},{MateChromosome},{MatebinPosition}"
    #write values to the bin key as defined in the line above
    Dict[BinKey] += [[A[0],int(A[1]),A[2],int(A[3])]]


    PreviousReadPosition = int(A[1])
    PreviousMatePosition = int(A[3])


with open(str(snakemake.output), 'w') as f:
    PrevDataSummary = ["chr",0,150,"chr",0,150]
    PrevKey = "Chr, 0, Chm, 0"
    for key, value in Dict.items():
        if len(value) >= 3:
            #Dict[key] = Elephant
            f.write('%s,%s,%s,%s' % (key, DataSummary(value)[0], value, Elephant))
            checkBin(DataSummary(value), key, PrevDataSummary,PrevKey, Elephant)
            PrevDataSummary = DataSummary(value)
            PrevKey = key
        else:
            pass
            #MAKE THIS PASS WHEN BUG CHECKING IS DONE
            #f.write('%s,%s,%s,%s' % (key, DataSummary(value), value, Elephant))

