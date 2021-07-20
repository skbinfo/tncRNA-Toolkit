#! /usr/bin/env python3
'''
This script do the job of pairwaise alignment
'''
import pandas as pd, sys
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

df1 = pd.read_csv(sys.argv[1], sep='\t')
df2 = pd.read_csv(sys.argv[2], sep='\t')
df3 = df1.merge(df2, on='tRNA_info')

for i in range(len(df3)):
    x=df3.loc[i, "tseq"]
    y=df3.loc[i, "Sequence"]
    for j,k in zip(x.split(),y.split()):
        alignments = pairwise2.align.localms(j, k, 1, -1, -1, -.5)
        for a in alignments:
            print('Type:'+df3.loc[i, 'tncRNA_Type'], '  Origin = '+df3.loc[i, 'tRNA_info'], '  Position:'+df3.loc[i, 'Position(start-end)']+'('+df3.loc[i,'Strand']+')', '  Length:'+str(df3.loc[i,'Length']))
            print(format_alignment(*a,full_sequences=True))
