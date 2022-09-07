import sys, os
'''
Adding zeroes to scaffold names
Written by Arielle Johnson 5/14/2022
Python 3
'''
print("INP1: FASTA file")
print("INP2: outfile name")
print("Output is same FASTA with 0s added to scaffold number. Assumes scaffold number after 2nd underscore.")

#open file and make out file
file1=open(sys.argv[1],'r')
out1=open(sys.argv[2],'w')
line1=file1.readline()
while line1:
    #print(line1)
    if line1.startswith(">"):
        split1=line1.strip().split('_')
        old_scaf_number=split1[2]

        scaf_number=str(old_scaf_number).zfill(4)

        new_scaf_name=split1[0]+"_"+split1[1]+"_"+scaf_number

        out1.write('{}\n'.format(new_scaf_name))

    else:
        out1.write(line1)

    line1=file1.readline()
file1.close()
out1.close()

print("Done")
