

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
import subprocess
from collections import defaultdict

#'python ' + scripts + 'SV_FindTEs.py {input.SV} {input.ref} {input.DB} {params.SV_Type} {wildcards.chromosome} {output.outfile}'

SVData = str(sys.argv[1])
REF = str(sys.argv[2])
DB = str(sys.argv[3])

SV = str(sys.argv[4])
CHROMin = str(sys.argv[5])

Output = str(sys.argv[6])  

TE_OVERLAP_THRESHOLD = 1000

dataout = []

print(CHROMin)

with open(Output, "w") as f:
    writer = csv.writer(f)

    if SV == "TandemDups":
        with open(SVData , "rt", newline = '') as file:
            rows = csv.reader(file, delimiter=',')
            data = list(rows)
            writer.writerow(data[0] + ["TE_hits","TE_s","TE_e"])
            #dataout.append(data[0] + ["TE_hits","TE_s","TE_e"])
            for row in data[1:]:
                Chrom = row[0]
                if CHROMin == Chrom:
                    Break_Points_Hit = 0
                    TE_names_s = "none"
                    TE_names_e = "none"
                

                    if row[5] != "Fail": 

                        subprocesstorun = "samtools faidx " + REF + " " + Chrom +':'+ str(int(row[3])-1000)+'-'+str(int(row[3])+1000) + "| blastn -db " + DB + " -max_target_seqs 1 -outfmt 6"

                        BlastOut = subprocess.check_output(subprocesstorun, shell=True,universal_newlines=True)
                        if len(BlastOut) > 0:
                            BlastOut = BlastOut.split()
                            if float(BlastOut[2]) > 80:
                                Break_Points_Hit += 1
                                TE_names_s = BlastOut[1]
                        
                        subprocesstorun = "samtools faidx " + REF + " " + Chrom +':'+ str(int(row[4])-1000)+'-'+str(int(row[4])+1000) + "| blastn -db " + DB + " -max_target_seqs 1 -outfmt 6"
                        BlastOut = subprocess.check_output(subprocesstorun, shell=True,universal_newlines=True)
                        if len(BlastOut) > 0:
                            BlastOut = BlastOut.split()
                            if float(BlastOut[2]) > 80:
                                Break_Points_Hit += 10
                                TE_names_e = BlastOut[1]
                
                                    

                    writer.writerow(row + [Break_Points_Hit,TE_names_s,TE_names_e])

    


    elif SV == "Rearrangements":
        with open(SVData , "rt", newline = '') as file:
            rows = csv.reader(file, delimiter=',')
            data = list(rows)
            writer.writerow(data[0] + ["TE_hits","TE_p","TE_m"])
            for row in data[1:]:
                Chrom = row[0]
                if CHROMin == Chrom:
                    Break_Points_Hit = 0
                    TE_names_p = "none"
                    TE_names_m = "none"
                    
                    
                    if row[2] == "=":
                        Chrom2 = Chrom
                    else:
                        Chrom2 = row[2]


                    if row[9] != "Fail" and CHROMin == Chrom:
                        subprocesstorun = "samtools faidx " + REF + " " + Chrom +':'+ str(int(row[5])-1000)+'-'+str(int(row[6])+1000) + "| blastn -db " + DB + " -max_target_seqs 1 -outfmt 6"
                        BlastOut = subprocess.check_output(subprocesstorun, shell=True,universal_newlines=True)
                        if len(BlastOut) > 0:
                            BlastOut = BlastOut.split()
                            if float(BlastOut[2]) > 80:
                                Break_Points_Hit += 1
                                TE_names_p = BlastOut[1]
                        
                        subprocesstorun = "samtools faidx " + REF + " " + Chrom2 +':'+ str(int(row[7])-1000)+'-'+str(int(row[8])+1000) + "| blastn -db " + DB + " -max_target_seqs 1 -outfmt 6"
                        BlastOut = subprocess.check_output(subprocesstorun, shell=True,universal_newlines=True)
                        if len(BlastOut) > 0:
                            BlastOut = BlastOut.split()
                            if float(BlastOut[2]) > 80:
                                Break_Points_Hit += 10
                                TE_names_m = BlastOut[1]
            
                                                    
                    writer.writerow(row + [Break_Points_Hit,TE_names_p,TE_names_m])

                    


        
     



