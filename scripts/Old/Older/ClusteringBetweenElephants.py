import re
import csv
import ast

OutputDic = {}

RangeDic = {}

INPUT = open(str(snakemake.input), 'r')
SV = str(snakemake.wildcards[0])

if SV == "Rearrangements":
    for line in INPUT:
        
        A = line.split(",")
        
        Info = ",".join(A[0:4])
            
        ElephantName = re.sub("\n", '', A[-1])

        Coverage = (len(A)-13)/4 #this is a lil weird - so the "line.split" code splits up the coverage list into a bunch of items. So we just find the number of items an entry would have with 0 coverage (len(A)-13) then divide it by the number of items in each coverage data (4)

        PairRANGE = [A[5], A[6]]
        MateRANGE = [A[9], A[10]]

        if Info in OutputDic:
            OutputDic[Info][1].append(ElephantName)
            OutputDic[Info][0].append(Coverage)


            RangeDic[Info][0] += PairRANGE
            RangeDic[Info][1] += MateRANGE

            MinPair = min(RangeDic[Info][0])
            MaxPair = max(RangeDic[Info][0])

            PairRange = [MinPair, MaxPair]

            MinMate = min(RangeDic[Info][1])
            MaxMate = max(RangeDic[Info][1])

            MateRange = [MinMate, MaxMate]

            RangeDic[Info] = [PairRange, MateRange]



        else:
            OutputDic[Info] = [[Coverage],[ElephantName]]
            RangeDic[Info] = [PairRANGE,MateRANGE] # Info = [["y","y"],["y","y"]]

    data = []

    for key, value in OutputDic.items():
        KeyValues = key.split(',')
        data.append([KeyValues[0],int(KeyValues[1]),KeyValues[2],int(KeyValues[3]), len(value[1]), RangeDic[key][0][0],RangeDic[key][0][1],RangeDic[key][1][0],RangeDic[key][1][1],value[0], value[1]])
        #f.write('%s,%s,%s,%s,%s' % (key, len(value), RangeDic[key], value, "\n"))

    startData = 0 
    endData = 1

    while startData != endData: 
        startData = data
        for i in range(len(data)):
            #checking if they're in the same chromosome
            try:
                if data[i][0] == data[i+1][0] and data[i][2] == data[i+1][2]:
                    #Checking for overlaps
                    Pair_Overlap = len(range(max(int(data[i][5]), int(data[i+1][5])), min(int(data[i][6]), int(data[i+1][6]))+1))
                    Mate_Overlap = len(range(max(int(data[i][7]), int(data[i+1][7])), min(int(data[i][8]), int(data[i+1][8]))+1))
                    if Pair_Overlap > 0 and Mate_Overlap > 0:
                        Overlap = True
                    else:
                        Overlap = False
                else:
                    Overlap = False
                    
                if Overlap == True:
                    #NewElephantNames =  [*set(data[i][-int(data[i][4]):] + data[i+1][-int(data[i+1][4]):])]
                    NewElephantNames = list(set(data[i][-1]+data[i+1][-1]))
                    NewCoverageList = data[i][-2]+data[i+1][-2]
                    pairRange = [min(data[i][5], data[i+1][5]) , max(data[i][6], data[i+1][6])]
                    pairRange.sort()

                    mateRange = [min(data[i][7], data[i+1][7])] + [max(data[i][8], data[i+1][8])]
                    mateRange.sort()

                    data[i+1] = data[i][0:4] + [len(NewElephantNames)] + pairRange + mateRange + [NewCoverageList] + [NewElephantNames]
                    data.remove(data[i])

                endData = data
            
            except IndexError:
                endData = data
                #print("end", len(data))
                break


    ###THIS IS A FUNCTION TO FILTER OUT ANY MUTATIONS THAT AREN'T FOUND AT LEAST 3 TIMES IN AT LEAST 1 ELEPHANT
    for i in range(len(data)):
        if max(data[i][-2]) >= 3:
            data[i][-2] = sum(data[i][-2])/len(data[i][-2])
        else:
            data[i][-2] = "Fail"

    header =  ["Chromosome","ReadBinPosition","MateChromosome","MateBinPosition","NumberOfElephants","PairMin","PairMax","MateMin","MateMax","CoverageList","ElephantNames"]
    data.insert(0,header)
    with open(str(snakemake.output), 'w', newline="") as f:
        writer = csv.writer(f)
        writer.writerows(data)

if SV == "TandemDups":
    for line in INPUT:
    
        A = line.split(",")
        
        Info = ",".join(A[0:2])
            
        ElephantName = re.sub("\n", '', A[-1])
        Coverage = (len(A)-7)/3 #this is a lil weird - so the "line.split" code splits up the coverage list into a bunch of items. So we just find the number of items an entry would have with 0 coverage (len(A)-7) then divide it by the number of items in each coverage data (3)
        PairRANGE = [A[3], A[4]]

        if Info in OutputDic:
            OutputDic[Info][1].append(ElephantName)
            OutputDic[Info][0].append(Coverage)


            RangeDic[Info][0] += PairRANGE

            MinPair = min(RangeDic[Info][0])
            MaxPair = max(RangeDic[Info][0])

            PairRange = [MinPair, MaxPair]

            RangeDic[Info] = [PairRange]



        else:
            OutputDic[Info] = [[Coverage],[ElephantName]]
            RangeDic[Info] = [PairRANGE]

    data = []

    for key, value in OutputDic.items():
        KeyValues = key.split(',')
        data.append([KeyValues[0],int(KeyValues[1]), len(value[1]), RangeDic[key][0][0],RangeDic[key][0][1],value[0], value[1]])
        #f.write('%s,%s,%s,%s,%s' % (key, len(value), RangeDic[key], value, "\n"))

    startData = 0 
    endData = 1

    while startData != endData: 
        startData = data
        for i in range(len(data)):
            try:
                #checking if they're in the same chromosome  
                if data[i][0] == data[i+1][0]:
                    #Checking for overlaps
                    Pair_Overlap = len(range(max(int(data[i][3]), int(data[i+1][3])), min(int(data[i][4]), int(data[i+1][4]))+1))
                    if Pair_Overlap > 0:
                        Overlap = True
                    else:
                        Overlap = False
                else:
                    Overlap = False
                    
                if Overlap == True:
                    #NewElephantNames =  [*set(data[i][-int(data[i][4]):] + data[i+1][-int(data[i+1][4]):])]
                    NewElephantNames = list(set(data[i][-1]+data[i+1][-1]))
                    NewCoverageList = data[i][-2]+data[i+1][-2]
                    pairRange = [min(data[i][3], data[i+1][3]) , max(data[i][4], data[i+1][4])]
                    pairRange.sort()

                    data[i+1] = data[i][0:2] + [len(NewElephantNames)] + pairRange + [NewCoverageList] + [NewElephantNames]
                    data.remove(data[i])

                endData = data
            
            except IndexError:
                endData = data
                #print("end", len(data))
                break


    ###THIS IS A FUNCTION TO FILTER OUT ANY MUTATIONS THAT AREN'T FOUND AT LEAST 3 TIMES IN AT LEAST 1 ELEPHANT
    for i in range(len(data)):
        if max(data[i][-2]) >= 3:
            data[i][-2] = sum(data[i][-2])/len(data[i][-2])
        else:
            data[i][-2] = "Fail"

    header =  ["Chromosome","ReadBinPosition","NumberOfElephants","PairMin","PairMax","CoverageList","ElephantNames"]
    data.insert(0,header)
    with open(str(snakemake.output), 'w', newline="") as f:
        writer = csv.writer(f)
        writer.writerows(data)