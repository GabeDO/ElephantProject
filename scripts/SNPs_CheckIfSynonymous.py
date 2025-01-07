from collections import defaultdict
import subprocess
import csv
import gzip


#setting up a dictionary for amino acids.
CodonDIC = defaultdict(list)
CodonDIC['I'] = ['ATT','ATC','ATA']
CodonDIC['L'] = ['CTT','CTC','CTA', 'CTG', 'TTA', 'TTG']
CodonDIC['V'] = ['GTT', 'GTC', 'GTA', 'GTG']
CodonDIC['F'] = ['TTT', 'TTC']
CodonDIC['M'] = ['ATG']
CodonDIC['C'] = ['TGT', 'TGC']
CodonDIC['A'] = ['GCT', 'GCC', 'GCA', 'GCG']
CodonDIC['G'] = ['GGT', 'GGC', 'GGA', 'GGG']
CodonDIC['P'] = ['CCT', 'CCC', 'CCA', 'CCG']
CodonDIC['T'] = ['ACT', 'ACC', 'ACA', 'ACG']
CodonDIC['S'] = ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC']
CodonDIC['Y'] = ['TAT', 'TAC']
CodonDIC['W'] = ['TGG']
CodonDIC['Q'] = ['CAA', 'CAG']
CodonDIC['N'] = ['AAT', 'AAC']
CodonDIC['H'] = ['CAT', 'CAC']
CodonDIC['E'] = ['GAA', 'GAG']
CodonDIC['D'] = ['GAT', 'GAC']
CodonDIC['K'] = ['AAA', 'AAG']
CodonDIC['R'] = ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG']
CodonDIC['Stop'] = ['TAA', 'TAG', 'TGA']

def CheckSNP(REF, SNP):
    #iterate through the Codon dictionary
    for key, value in CodonDIC.items():
        RefCheck = ''
        SNPCheck = ''
        #if the reference or SNP'd codon is in the amino acid group, flag it 
        if REF in value:
            RefCheck = key
        if SNP in value:
            SNPCheck = key
        #if either of the things are flaged, do some shit
        if RefCheck != '' or SNPCheck != '':
            #check for a missing value, and if its there, see if the two values can possilbly mathc
            if 'N' in SNP:
                #if they can't match no matter what, its NS
                if bool(set(CanItFit(SNP)) & set(value)) == False:
                    return 'NS'
                else:
                    return 'Sn*'
            #same thing for if theres a missing value in the REF codon
            if 'N' in REF:
                if bool(set(CanItFit(REF)) & set(value)) == False:
                    return 'NS'
                else:
                    return 'Sn*'
            
            #if theres no missing values, see if bioth the SNP and the REF have been found in that codon group - if not, they're NS
            if RefCheck == SNPCheck:
                return 'Sn'
            else:
                return 'NS'
            
        #if they both have a missibg value, and are therefore not found in the codon dict  
        if 'N' in SNP:
            if bool(set(CanItFit(SNP)) & set(value)) == False:
                return 'NS'
            else:
                return 'Sn*'
        
        if 'N' in REF:
            if bool(set(CanItFit(REF)) & set(value)) == False:
                return 'NS'
            else:
                return 'Sn*'
    
    
def CanItFit(Seq):
    #set up some lists for later
    PermList = []
    UnKnownPos =[]
    NList = ['A','T','C','G']

    #loop through the sequence positions
    for i in range(len(Seq)):
        #if theres a missing value, append that position to a list of positions
        if Seq[i] == 'N':
            UnKnownPos.append(i)
    
    #if theres only a single missing value, just loop through ass possibilities and add em to a list
    if len(UnKnownPos) == 1:
        for i in range(len(NList)):
            Seq1 = Seq[:UnKnownPos[0]] + NList[i] + Seq[UnKnownPos[0] + 1:]
            PermList.append(Seq1)
        
    #same thing for if theres two, but nested loop it, so it gets all possible permutaions. 
    if len(UnKnownPos) == 2:
        for i in range(len(NList)):
            Seq1 = Seq[:UnKnownPos[0]] + NList[i] + Seq[UnKnownPos[0] + 1:]
            for k in range(len(NList)):
                Seq2 = Seq1[:UnKnownPos[1]] + NList[k] + Seq1[UnKnownPos[1] + 1:]
                
                PermList.append(Seq2)

    #return a list of all possibilities       
    return PermList
        

def GetCodon(pos, Sequence,snp):
    CorSeq = list(Sequence)
    try:
        CorSeq.remove(' ')
    except ValueError:
        pass
    Codon = CorSeq[(pos - pos % 3):(pos - pos % 3)+3]
    
    if snp:
        SNPpos_ = pos%3
        Codon[SNPpos_] = snp
    
    return ' '.join(Codon)




print("Making gene dictionary")
GTF = str(snakemake.input[0])
with gzip.open(GTF, "rt", newline = '') as file: ##with gzip.open (GTF, "rt", newline = ") as file:
    rows = csv.reader (file, delimiter='\t')
    for row in rows:
        if len(row) != 1 and row[8].split("\"")[5][0] != "(":
            if row[2] in ["gene", "CDS", "start_codon","stop_codon"]: #CHECK IF ITS ONE OF THE TYPES I CARE ABOUT
                ## I have to use 'row [0].replace(".", "_")' to turn the "." in
                ##the chrom names to "_" because python dosen't like when variables have "." in their names
                if row[0].replace(".", "_") in locals (): # IF THE CHROMDICTONARY EXISTS
                    #ADD TO THE DICTONARY IN THE FORMAT: CHROM= {GENEID: START, END, TYPE} 
                    eval(row[0].replace(".", "_"))[row[8].split("\"")[5]] += [ [int(row[3]),int(row[4]),row[2]] ]
                else: # IF THE CHROMDICTONARY DOSEN'T EXISTS
                    print(row[0].replace(".", "_")) # just prints out the Chrom name. If all is working, each chrom is printed once
                    exec(row[0].replace(".", "_") + "= defaultdict(list)") #MAKE A DICTONARY NAMED AFTER THE CHROMOSOME 
                    #ADD TO THE DICTONARY IN THE FORMAT: CHROM= {GENEID: START, END, TYPE}
                    eval(row[0].replace(".", "_"))[row[8].split("\"")[5]] += [ [int(row[3]),int(row[4]),row[2]] ]
print("Gene dictionary complete!")



"""
print("Making SNP dictionary")
SNPList = []
VCFfile = str(snakemake.input[1])
with gzip.open(VCFfile, 'rt') as file:
    rows = csv.reader(file, delimiter='\t')
    for row in rows:
        if row[0][0] != "#":
                AC = sum([ int(x) for x in (row[7]+"AC=0").split("AC=")[1].split()[0].split(",") ])
                AN = int((row[7]+"AN=0").split("AN=")[1].split()[0].split(";")[0])
                SNPList.append([row[0], int(row[1]), row[3], row[4], AC/AN ])##This just makes a list of all the SVsI have. The row 7 stuff is the Allele freq
print("SNP dictionary complete!")
"""

print("Checking SNPs against Genes.\nThis may take some time...")
ResultOutput = []
GENOfile = str(snakemake.input[1])
with open(GENOfile, 'rt') as file:
    rows = csv.reader(file, delimiter='\t')
    for i in rows:
        if i[0] != "CHROM":
            #get some info from the VCF file per row
            Chrom = i[0]
            SNPpos = int(i[1])
            REF = i[2]
            SNP_ = str(i[3])
            AlFreq = int(i[6]) * 2 + int(i[8])
            Polarized = i[9]

            GeneFound = False
            GeneFoundAdj = False
            
            Result = "NA"
            #check if the chrom has a GTF dict
            if Chrom.replace(".", "_") in locals ():
                #if it does, run through it untill it finds a CDS with the correct position.
                for key, values in eval(Chrom.replace(".", "_")).items():
                    #give it a defult result that can change later
                    CDS_Start_FIRST = None
                    CDS_End_FIRST = None
                    Gene_Start = -1
                    Gene_End = -1
                    TESTHIT = 0
                    for s in values:
                        if GeneFound == True:
                            break
                        if s[2] == "gene":
                            Gene_Start = s[0]
                            Gene_End = s[1]
                            if SNPpos >= s[0] and SNPpos <= s[1]:
                                Result = "On Gene not on CDS"
                            elif SNPpos >= s[0]-10000 and SNPpos <= s[1]+10000:
                                Result = "Adjacent To Gene"
                            elif Result == "NA":
                                Result = "Not Near Gene"
                        else:
                            if s[2] == "CDS":
                                if CDS_Start_FIRST == None:
                                    CDS_Start_FIRST = s[0]
                                if CDS_End_FIRST == None:
                                    CDS_End_FIRST = s[1]

                                if s[0] < CDS_Start_FIRST:
                                    CDS_Start_FIRST = s[0]
                                if s[1] > CDS_End_FIRST:
                                    CDS_End_FIRST = s[1]
                                
                                
                            #If the SNP position is found at one of the CDS sites
                            if SNPpos >= s[0] and SNPpos <= s[1]:
                                GeneFound = True
                                CDS_Start = s[0]
                                CDS_End = s[1]
                                #Set up a sub process, this will use samtools to pull out the sequence as a string from the reference file
                                ProcessToRun_CDS_Seq = "echo $(samtools faidx " + str(snakemake.input[2]) + " " + str(Chrom) + ":" + str(CDS_Start) + "-" + str(CDS_End) + "| grep -v '>')"
                                Seq = subprocess.check_output(ProcessToRun_CDS_Seq, shell=True)
                                #the subprocess returns the sequence as bytes, but we want it in plaintext, so we gotta decode it (.decode('ascii')). The .strip() command removes the \n thats atted at the end of it. 
                                Seq = (Seq.decode('ascii')).strip()

                                #get the relative position of the SNP on the CDS site - the code in 'GetCodon' counts from 0, but the VCF/GFT file positions don't.
                                #I would normally have to add a +1 to the 'SNPpos - s' equation to get the exact position, but since pythoin counts from 0, the correct
                                #pos value will be 1 less than what the VCF/GTF say, so I just omit the +1 (I think thats correct anyway)
                                RelativeSNPpos = SNPpos - s[0]

                                #run the code
                                #for some reason, it was outputting the Codons a little weird, with a space between each neucliotide and a missing value as an extra space.
                                #So I'm just replacing 3 consecutive spaces with an 'N' and then removing all the other spaces to give it a nice clean codon
                                Ref = GetCodon(RelativeSNPpos,Seq,False)
                                Ref = Ref.replace("   ", "N")
                                Ref = Ref.replace(" ", "")

                                SNP = GetCodon(RelativeSNPpos,Seq,SNP_)
                                SNP = SNP.replace("   ", "N")
                                SNP = SNP.replace(" ", "")

                                
                                #here I'm just chcecking if the SNPS are synonomous for each possible neucliotide. This is to check if they're 4-fold Sn (eg: will be Sn nomatter what the mutation is)
                                #I append each result to the list 'ffns' and then later check if there's 0 NS values in there (eg, none of the possibilities are NS)
                                ffSn = []
                                for i in ['A','T','C','G']:
                                    FourFoldTest_iter = GetCodon(RelativeSNPpos,Seq,i)
                                    FourFoldTest_iter = FourFoldTest_iter.replace("   ", "N")
                                    FourFoldTest_iter = FourFoldTest_iter.replace(" ", "")
                                    ffSn.append(CheckSNP(Ref,FourFoldTest_iter))
                                

                                #write results to the output list, then break the loop so they code stopps looking (and dosen't trigger the else statements)
                                #this is where I check if its 4-fold Sn "ffns.count('NS') == 0", and if it is, I just tag is 4fSn. Otherwise, I'll check as normal. 
                                if ffSn.count('NS') == 0:
                                    Result = "Sn_4f"
                                elif ffSn.count('Sn')+ffSn.count('Sn*') == 0:
                                    Result = "NS_4f"
                                else:
                                    Result = CheckSNP(Ref,SNP)
                                break

                            elif SNPpos >= s[0]-10000 and SNPpos <= s[1]+10000 and GeneFound == False and Result != "UTR" and s[2] == "CDS":
                                Result = "Adjacent to CDS"
                                GeneFoundAdj = True

                            elif Result == "NA":
                                Result = "Not Near Gene"
                    #############################################
                    if CDS_Start_FIRST != None:
                        if SNPpos > Gene_Start and SNPpos < CDS_Start_FIRST:
                            Result = "UTR_s"

                    if CDS_End_FIRST != None:
                        if SNPpos < Gene_End and SNPpos > CDS_End_FIRST:
                            Result = "UTR_e"

                    #if the result is 
                    if GeneFound == True:
                        break
            
            if Result != "NA":
                ResultOutput.append([Chrom, SNPpos, Result,AlFreq,Polarized])

print("process complete!\nwritting results to file...")

Outputfile = str(snakemake.output[0])
with open(Outputfile, 'w') as f: 
    f.write('%s,%s,%s,%s,%s\n' % ("Chrom","Pos","Result","AlFreq","Pol"))
    #f.write('Chrom,Position,Result\n')
    for i in ResultOutput:
        if len(str(i[3]).split(",")) == 1:
            f.write('%s,%s,%s,%s,%s\n' % (i[0],i[1],i[2],i[3],i[4]))


Outputfile2 = str(snakemake.output[1])
with open(Outputfile2, 'w') as f: 
    f.write('%s,%s,%s,%s,%s\n' % ("Chrom","Pos","Result","AlFreq","Pol"))
    #f.write('Chrom,Position,Result\n')
    for i in ResultOutput:
        if i[2][0:2] == "NS" or i[2][0:2] == "Sn" and len(str(i[3]).split(",")) == 1:
            f.write('%s,%s,%s,%s,%s\n' % (i[0],i[1],i[2],i[3],i[4]))

print("All complete!")