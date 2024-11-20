from collections import defaultdict
import numpy as np
import csv
import sys
import gzip


data =[]
CovData = defaultdict(list)
VarData = defaultdict(list)
BinSize = 10000

CovData_mean = defaultdict(float)
VarData_mean = defaultdict(float)

#smINPUT = str(snakemake.input)

with gzip.open('CoverageAndVariance.out.gz', "rt") as lines:
        for i in lines:
            chom, pos, cov, var = i.split(" ")
            CovData[str(chom),int(pos)-int(pos)%BinSize] += [float(cov)]
            if len(CovData[str(chom),int(pos)-int(pos)%BinSize]) >= BinSize:
                CovData_mean[str(chom),int(pos)-int(pos)%BinSize] = np.mean(CovData[str(chom),int(pos)-int(pos)%BinSize])
                del CovData[str(chom),int(pos)-int(pos)%BinSize]

            VarData[str(chom),int(pos)-int(pos)%BinSize] += [float(var)]
            if len(VarData[str(chom),int(pos)-int(pos)%BinSize]) >= BinSize:
                VarData_mean[str(chom),int(pos)-int(pos)%BinSize] = np.mean(VarData[str(chom),int(pos)-int(pos)%BinSize])
                del VarData[str(chom),int(pos)-int(pos)%BinSize]


for i in CovData:
    CovData_mean[i] = np.mean(CovData[i])
for i in VarData:
    VarData_mean[i] = np.mean(VarData[i])

del CovData
del VarData


header = ['Chromosome','Position','Coverage','Variance']
with open('BinnedCoverage.csv', 'w') as f:
    write = csv.writer(f)
    write.writerow(header)
    for key, value in CovData_mean.items():
        
        Chromosome = key[0]
        Position = key[1]            
        
        write.writerow([Chromosome, Position, CovData_mean[key], VarData_mean[key]])
