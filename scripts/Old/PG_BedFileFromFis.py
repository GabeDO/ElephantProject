import csv

FISFILE = str(snakemake.input[0])
BinSize = 100000

BEDFILE_output_list = []

with open(FISFILE, newline='') as file:
    rows = csv.reader(file, delimiter=',')
    for row in rows:
        try:
            if float(row[2]) > 0.2:
                BEDFILE_output_list.append([row[0],int(row[1]),int(row[1])+BinSize-1])
        except ValueError:
            print('the header is: ' + str(row))


BEDFILE_output = str(snakemake.output[0])
with open(BEDFILE_output, 'w') as BEDFILE_output:
    wr = csv.writer(BEDFILE_output, delimiter='\t')
    wr.writerows(BEDFILE_output_list)