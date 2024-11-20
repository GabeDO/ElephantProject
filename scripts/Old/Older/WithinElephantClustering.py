from datetime import date
today = date.today()
import gzip
import sys
import re
import statistics
import os
from collections import defaultdict





def FixBin(Dictonary, key1, key2):
    Dictonary[key2] += Dictonary[key1]
    Dictonary.pop(key1)

SV=str(snakemake.wildcards[1])

if SV == 'Rearrangements':

    #function for checking if binning works, add " _"+str(today)+" " to the output file name if you wanna tag it with dates
    def checkBin(DataSum, key, PrevDataSum, prevkey, Elephant, output):
        with open("BinErrors.txt", 'a') as Err:  
            ReadError = ""
            MateError = ""
            if DataSum[0][0] == PrevDataSum[0][0] and DataSum[0][4] == PrevDataSum[0][4]:
                if PrevDataSum[0][1]-350 <= DataSum[0][1] <= PrevDataSum[0][2]+350 and PrevDataSum[0][5]-350 <= DataSum[0][5] <= PrevDataSum[0][6]+350:
                    if PrevDataSum[0][1] <= DataSum[0][1] <= DataSum[0][2] <= PrevDataSum[0][2]:
                        ReadError = "Read within previous"
                        ERROR = False
                    if PrevDataSum[0][5] <= DataSum[0][5] <= DataSum[0][6] <= PrevDataSum[0][6]:
                        MateError = "Mate within previous"
                        ERROR = False
                    if len(ReadError) + len(MateError) == 0:
                        if output == True:
                            Err.write(f"{key}:{DataSum[0]},Reads in bin:,{DataSum[1]} ||||| {prevkey}:{PrevDataSum[0]},prev reads in bin,{PrevDataSum[1]},{Elephant}")
                        ERROR = True
                    else:
                        ERROR = False
                else:
                    ERROR = False
            else:
                ERROR = False
        
        return [ERROR, [key,prevkey]]

    LoopCheckValue = 0
    Dict = defaultdict(list)
    binsize = 700
    ReadLen = 151
    OutputList = []

    PreviousReadPosition = 0
    PreviousMatePosition = 0
    INPUT = open(str(snakemake.input), 'r')

    PreviousBin = 0
    PreviousMateBin = 0

    PrevDataSummary = ["chr",0,150,"chr",0,150]
    PrevKey = "Chr, 0, Chm, 0"

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

    for line in INPUT:
        
        A=line.split(' ')
        Chromosome = A[0]
        ReadPosition = int(A[1])
        MateChromosome = A[2]
        MatePosition = int(A[3])
        Elephant = A[-1]
        
        ReadbinPosition = (ReadPosition-350) - (ReadPosition-350) % (binsize)
        PreviousBin = (ReadPosition-350) - (ReadPosition-350) % (binsize)

        PreviousMateBin = (MatePosition-350) - (MatePosition-350) % (binsize)
        MatebinPosition = (MatePosition-350) - (MatePosition-350) % (binsize)
        
        BinKey=f"{Chromosome},{ReadbinPosition},{MateChromosome},{MatebinPosition}"
        Dict[BinKey] += [[A[0],int(A[1]),A[2],int(A[3])]]


        PreviousReadPosition = int(A[1])
        PreviousMatePosition = int(A[3])
        PreviousBin = (ReadPosition-350) - (ReadPosition-350) % (binsize)
        PreviousMateBin = (MatePosition-350) - (MatePosition-350) % (binsize)


    newDict = "yeet"
    while newDict != list(Dict):
        newDict = list(Dict)
        LoopCheckValue += 1
        for key in list(Dict):
            DS = DataSummary(Dict[key])
            x = checkBin(DS, key, PrevDataSummary,PrevKey, Elephant, False)
            if x[0] == True:
                FixBin(Dict, x[1][0], x[1][1])    
            PrevDataSummary = DS
            PrevKey = key


    with open(str(snakemake.output), 'w') as f:
        PrevDataSummary = ["chr",0,150,"chr",0,150]
        PrevKey = "Chr, 0, Chm, 0"
        for key, value in Dict.items():
            if len(value) >= 2:
                #Dict[key] = Elephant
                f.write('%s,%s,%s,%s' % (key, DataSummary(value)[0], value, Elephant))
                checkBin(DataSummary(value), key, PrevDataSummary,PrevKey, Elephant, True)
                PrevDataSummary = DataSummary(value)
                PrevKey = key
            else:
                checkBin(DataSummary(value), key, PrevDataSummary,PrevKey, Elephant, True)
                PrevDataSummary = DataSummary(value)
                PrevKey = key
                #MAKE THIS PASS WHEN BUG CHECKING IS DONE
                #f.write('%s,%s,%s,%s' % (key, DataSummary(value), value, Elephant))


if SV == 'TandemDups':

    #function for checking if binning works, add " _"+str(today)+" " to the output file name if you wanna tag it with dates
    def checkBin(DataSum, key, PrevDataSum, prevkey, Elephant, output):
        with open("BinErrors.txt", 'a') as Err:  
            ReadError = ""
            MateError = ""
            if DataSum[0][0] == PrevDataSum[0][0]:
                if PrevDataSum[0][1]-350 <= DataSum[0][1] <= PrevDataSum[0][2]+350:
                    if PrevDataSum[0][1] <= DataSum[0][1] <= DataSum[0][2] <= PrevDataSum[0][2]:
                        ReadError = "Read within previous"
                        ERROR = False
                    if len(ReadError) + len(MateError) == 0:
                        if output == True:
                            Err.write(f"{key}:{DataSum[0]},Reads in bin:,{DataSum[1]} ||||| {prevkey}:{PrevDataSum[0]},prev reads in bin,{PrevDataSum[1]},{Elephant}")
                        ERROR = True
                    else:
                        ERROR = False
                else:
                    ERROR = False
            else:
                ERROR = False
        
        return [ERROR, [key,prevkey]]

    LoopCheckValue = 0
    Dict = defaultdict(list)
    binsize = 700
    ReadLen = 151
    OutputList = []

    PreviousReadPosition = 0
    PreviousMatePosition = 0
    INPUT = open(str(snakemake.input), 'r')

    PreviousBin = 0
    PreviousMateBin = 0

    PrevDataSummary = ["chr",0,150,150,[0,150]]
    PrevKey = "Chr, 0"

    def DataSummary(LIST):
        ListTD = []
        TDstart = min([LIST[0][1],LIST[0][2]])
        TDend = max([LIST[0][1],LIST[0][2]])+ReadLen

        for i in LIST:
            ListTD.append(i[1])
            if (i[1] or i[2]) < TDstart:
                TDstart = min([i[1],i[2]])
            elif (i[1]+ReadLen or i[2]+ReadLen) > TDend:
                TDend = max([i[1],i[2]])+ReadLen
            else:
                pass

        return [[LIST[0][0],TDstart,TDend,(TDend-TDstart)],[ListTD]]


    for line in INPUT:
        
        A=line.split(' ')
        Chromosome = A[0]
        PositionStart = min([int(A[1]),int(A[3])])
        PositionEnd = min([int(A[1]),int(A[3])])
        Elephant = A[-1]
        
        binPosition = (PositionStart-350) - (PositionStart-350) % (binsize)

        BinKey=f"{Chromosome},{binPosition}"
        Dict[BinKey] += [[Chromosome,PositionStart,PositionEnd]]

        PreviousPositionStart = min([int(A[1]),int(A[3])])
        PreviousBin = (PositionStart-350) - (PositionStart-350) % (binsize)

    newDict = "yeet"
    while newDict != list(Dict):
        newDict = list(Dict)
        LoopCheckValue += 1
        for key in list(Dict):
            DS = DataSummary(Dict[key])
            x = checkBin(DS, key, PrevDataSummary,PrevKey, Elephant, False)
            if x[0] == True:
                FixBin(Dict, x[1][0], x[1][1])    
            PrevDataSummary = DS
            PrevKey = key


    with open(str(snakemake.output), 'w') as f:
        PrevDataSummary = ["chr",0,150,"chr",0,150]
        PrevKey = "Chr, 0, Chm, 0"
        for key, value in Dict.items():
            if len(value) >= 2:
                #Dict[key] = Elephant
                f.write('%s,%s,%s,%s' % (key, DataSummary(value)[0], value, Elephant))
                checkBin(DataSummary(value), key, PrevDataSummary,PrevKey, Elephant, True)
                PrevDataSummary = DataSummary(value)
                PrevKey = key
            else:
                checkBin(DataSummary(value), key, PrevDataSummary,PrevKey, Elephant, True)
                PrevDataSummary = DataSummary(value)
                PrevKey = key
                #MAKE THIS PASS WHEN BUG CHECKING IS DONE
                #f.write('%s,%s,%s,%s' % (key, DataSummary(value), value, Elephant))
