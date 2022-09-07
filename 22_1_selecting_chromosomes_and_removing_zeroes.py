import sys, os
'''
Selecting the 8 Euphorbia peplus chromosomes in a fasta file
Written by Arielle Johnson 5/18/22
Python 3
'''

#open file and make out file
file1=open(sys.argv[1],'r')
out1=open(sys.argv[1]+"_chromosomes.fasta",'w')

do_write=0
line1=file1.readline()
while line1:
    if line1.startswith(">"):
        #print(line1)
        split1=line1.strip().split('_')
        old_scaf_number=split1[1].replace("sc","")
        new_scaf_number=old_scaf_number.lstrip("0")
        new_scaf_name=split1[0]+"_chr"+new_scaf_number+"_"+split1[2]
        if int(new_scaf_number)<=8:
            out1.write('{}\n'.format(new_scaf_name))
            do_write=1

        else:
            do_write=0

    else:
        if do_write==1:
            out1.write(line1)
        elif do_write==0:
            pass

    line1=file1.readline()
file1.close()
out1.close()

print("Done")
