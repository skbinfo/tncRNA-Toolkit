#! /usr/bin/env python3
'''
This script count the number for each tncRNA class
'''
import pandas as pd, sys, os

file1=sys.argv[1]
if os.path.isfile(file1) and os.path.getsize(file1) != 0:
    print('tncRNA_Type:\tNumber')
    df = pd.read_csv(file1, sep='\t')
## Counting tsRNAs
    tRF_5 = df[df['tncRNA_Type']=="tRF-5"]
    if tRF_5.empty == False:
        tRF_5_count = pd.Series(tRF_5['tncRNA_Type']).count()
        print('tRF-5:\t'+str(tRF_5_count))
    else:
        print('tRF-5:\t0')

    half_5 = df[df['tncRNA_Type']=="5'tRH"]
    if half_5.empty == False:
        half_5_count = pd.Series(half_5['tncRNA_Type']).count()
        print('5\'tRH:\t'+str(half_5_count))
    else:
        print('5\'tRH:\t0')

    tRF_3CCA = df[df['tncRNA_Type']=="tRF-3-CCA"]
    if tRF_3CCA.empty == False:
        tRF_3CCA_count = pd.Series(tRF_3CCA['tncRNA_Type']).count()
        print('tRF-3-CCA:\t'+str(tRF_3CCA_count))
    else:
        print('tRF-3-CCA:\t0')
    
    tRF_3 = df[df['tncRNA_Type']=="tRF-3"]
    if tRF_3.empty == False:
        tRF_3_count = pd.Series(tRF_3['tncRNA_Type']).count()
        print('tRF-3:\t'+str(tRF_3_count))
    else:
        print('tRF-3:\t0')
    
    half_3CCA = df[df['tncRNA_Type']=="3'tRH-CCA"]
    if half_3CCA.empty == False:
        half_3CCA_count = pd.Series(half_3CCA['tncRNA_Type']).count()
        print('3\'tRH-CCA: '+str(half_3CCA_count))
    else:
        print('3\'tRH-CCA:\t0')


    half_3 = df[df['tncRNA_Type']=="3'tRH"]
    if half_3.empty == False:
        half_3_count = pd.Series(half_3['tncRNA_Type']).count()
        print('3\'tRH: '+str(half_3_count))
    else:
        print('3\'tRH:\t0')

    tRF_1 = df[df['tncRNA_Type']=="tRF-1"]
    if tRF_1.empty == False:
        tRF_1_count = pd.Series(tRF_1['tncRNA_Type']).count()
        print('tRF-1:\t'+str(tRF_1_count))
    else:
        print('tRF-1:\t0')
    
    ld_tRF = df[df['tncRNA_Type']=="leader-tRF"]
    if ld_tRF.empty == False:
        ld_tRF_count = pd.Series(ld_tRF['tncRNA_Type']).count()
        print('leader-tRF:\t'+str(ld_tRF_count))
    else:
        print('leader-tRF:\t0')

    i_tRF = df[df['tncRNA_Type']=="other-tRF"]
    if i_tRF.empty == False:
        i_tRF_count = pd.Series(i_tRF['tncRNA_Type']).count()
        print('other-tRF:\t'+str(i_tRF_count))
    else:
        print('other-tRF:\t0')

else:
    pass
