import sys, os
'''
Removing the gene and transcript feature lines in a gtf
Written by Arielle Johnson 12/17/2021
Python 3
'''
print("INP1: .gtf file")
print("Output is same gtf with gene and transcript feature lines removed")

#open file and make out file
file1=open(sys.argv[1],'r')
out1=open(sys.argv[1].replace('.gtf','')+"_no_gene_or_transcript.gtf",'w')
line1=file1.readline()
while line1:
    #print(line1)
    tab1=line1.strip().split('\t')
    feature=tab1[2]

    if feature=="gene":
        pass

    elif feature=="transcript":
        pass

    else:
        out1.write(line1)

    line1=file1.readline()
file1.close()
out1.close()

print("Done")
