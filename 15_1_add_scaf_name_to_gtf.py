import sys, os
'''
Adding scaffold names to a gtf (used TSEBRA output that was then renamed using rename_gtf.py)
Written by Arielle Johnson 12/10/2021, modified 5/11/22 to add zeroes
Python 3
'''
print("INP1: .gtf file from TSEBRA")
print("Output is same gtf with gene and transcript names changed to include scaffold number. Also adds zeroes to scaffold number and chromosome number.")

#open file and make out file
file1=open(sys.argv[1],'r')
out1=open(sys.argv[1].replace('.gtf','')+"_withScaf.gtf",'w')
line1=file1.readline()
while line1:
    #print(line1)
    tab1=line1.strip().split('\t')
    seqname=tab1[0]
    source=tab1[1]
    feature=tab1[2]
    start=tab1[3]
    end=tab1[4]
    score=tab1[5]
    strand=tab1[6]
    frame=tab1[7]
    attributes=tab1[8]

    old_scaf_number=seqname.strip().split('_')[2]

    if len(old_scaf_number)==4:
        scaf_number=old_scaf_number

    elif len(old_scaf_number)==3:
        scaf_number="0"+old_scaf_number

    elif len(old_scaf_number)==2:
        scaf_number="00"+old_scaf_number

    elif len(old_scaf_number)==1:
        scaf_number="000"+old_scaf_number

    new_seqname="HiC_scaffold_"+scaf_number

    if feature=="gene":
        gene_name_strip=attributes.strip().split("_")
        gene_and_scaf=gene_name_strip[0]+"_sc"+scaf_number
        new_attributes="{}_{}".format(gene_and_scaf, gene_name_strip[1])

    elif feature=="transcript":
        transcript_name_strip=attributes.strip().split("_")
        transcript_and_scaf=transcript_name_strip[0]+"_sc"+scaf_number
        new_attributes="{}_{}".format(transcript_and_scaf, transcript_name_strip[1])
    else:
        stripped_attributes=attributes.strip().split("\"")
        transcript_id=stripped_attributes[1]
        transcript_name_strip=transcript_id.strip().split("_")
        transcript_and_scaf=transcript_name_strip[0]+"_sc"+scaf_number
        old_transcript=transcript_name_strip[1].strip().split(".")
        old_transcript_number=old_transcript[0].replace('g','')

        if len(old_transcript_number)==5:
            new_transcript_number=old_transcript_number

        elif len(old_transcript_number)==4:
            new_transcript_number="0"+old_transcript_number

        elif len(old_transcript_number)==3:
            new_transcript_number="00"+old_transcript_number

        elif len(old_transcript_number)==2:
            new_transcript_number="000"+old_transcript_number

        elif len(old_transcript_number)==1:
            new_transcript_number="0000"+old_transcript_number

        new_transcript_name="g"+new_transcript_number+"."+old_transcript[1]

        new_transcript="{}_{}".format(transcript_and_scaf, new_transcript_name)

        gene_id=stripped_attributes[3]
        gene_name_strip=gene_id.strip().split("_")
        gene_and_scaf=gene_name_strip[0]+"_sc"+scaf_number
        old_gene_number=gene_name_strip[1].replace('g','')

        if len(old_gene_number)==5:
            new_gene_number=old_gene_number

        elif len(old_gene_number)==4:
            new_gene_number="0"+old_gene_number

        elif len(old_gene_number)==3:
            new_gene_number="00"+old_gene_number

        elif len(old_gene_number)==2:
            new_gene_number="000"+old_gene_number

        elif len(old_gene_number)==1:
            new_gene_number="0000"+old_gene_number

        new_gene_name="g"+new_gene_number

        new_gene="{}_{}".format(gene_and_scaf, new_gene_name)
        new_attributes="{}\"{}\"{}\"{}\";".format(stripped_attributes[0], new_transcript, stripped_attributes[2], new_gene)

    out1.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(new_seqname, source, feature, start, end, score, strand, frame, new_attributes))

    line1=file1.readline()
file1.close()
out1.close()

print("Done")
