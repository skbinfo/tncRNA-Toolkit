#! /usr/bin/env python3
'''
This script classify the tncRNA based on position and length
'''
import pandas as pd, sys, os
from datetime import datetime

# tRF-3-5-other-half
file1 = sys.argv[1] #bowtie output
file2a = file1+'.1mt.tmp'
file2b = file1+'.1pre.tmp'
file3 = file1+'.2tmp'
file4 = file1+'.3tmp'
s_cut = int(sys.argv[2]) #cutoff read-count

if os.path.isfile(file1) and os.path.getsize(file1) != 0:
    data = pd.read_csv(file1, header=None, sep='\t') #read file without header info and sep tab
    data.columns = ["qid", "strand", "sid", "start", "seq", "qual", "alignment", "mismatch"] #add header
    data.drop(["qual","alignment","mismatch"], axis = 1, inplace = True) #drop extra col
    data["sid"]= data["sid"].str.replace("\\(-\\)", "(#)")
    data["sid"]= data["sid"].str.replace("::|:", "\t")
    nu1 = data["start"]+data["seq"].str.len() #add start position with seq lenght to get end position
    data["start"]= data["start"]+1 #add +1 to lefter join offet 0
    data.insert(4, "end", nu1, True) #add nu1 to dataframe as end position
    nu2 = data["seq"].str.len() #get length of sequnce in nu2
    data.insert(6, "seq_len", nu2, True) #add seq_len
    
    dpre = data.loc[(data['sid'].str.contains('^ld_|^tl_'))] #separate dataframe for tRNA flank 50
    dmt = data.loc[(data['sid'].str.contains('^tRNA'))] #separate dataframe for mature tRNA

    dmt.to_csv(file2a, header= False, index= False, quoting=3, sep = '\t', escapechar=' ') #write dataframe in file
    dpre.to_csv(file2b, header= False, index= False, quoting=3, sep = '\t', escapechar=' ')

    dt = pd.read_csv(file2a, header=None, sep='\t')
    dt[4]= dt[4].str.replace("-|\\(", "\t")
    dt[4]= dt[4].str.replace("\\)", "")
    dt[4]= dt[4].str.replace("#", "-")
    dt.to_csv(file3, header= False, index= False, quoting=3, sep = '\t', escapechar=' ')
    
    df = pd.read_csv(file3, header=None, sep='\t')
    #nuu1 = df[5]-df[4]+1
    nuu1 = df[5]-df[4]
    df.insert(11, "tRNA_len", nuu1, True)

    #create reverse complementary if strand is -ve
    com_ = df[9].str.translate(str.maketrans('ATCG','TAGC'))
    rev_com = com_.str[::-1]
    df.loc[df[1] == '-', 9] = rev_com

    pos = df['tRNA_len']+3
    df.loc[(df[7]==1) & (df[10]>=14) & (df[10]<=30), 'category'] = 'tRF-5'
    df.loc[(df[7]==1) & (df[10]>30) & (df[10]<=50), 'category'] = "5'tRH"

    df.loc[(df[8]==pos) & (df[10]>=14) & (df[10]<=30) , 'category'] = 'tRF-3-CCA'
    df.loc[(df[8]==df['tRNA_len']) & (df[10]>=14) & (df[10]<=30) , 'category'] = 'tRF-3'

    df.loc[(df[8]==pos) & (df[10]>30) & (df[10]<=50), 'category'] = "3'tRH-CCA"
    df.loc[(df[8]==df['tRNA_len']) & (df[10]>30) & (df[10]<=50), 'category'] = "3'tRH"
    
    df.loc[(df[7]>1) & (df[8]<df['tRNA_len']) & (df[10]>=14) & (df[10]<=50), 'category'] = 'other-tRF'
    df.loc[(df[7]>1) & (df[8]>df['tRNA_len']) & (df[8]<pos) & (df[10]>=14) & (df[10]<=50), 'category'] = 'other-tRF'

    tloc = df[2]+'::'+df[3]+':'+df[4].astype(str)+'-'+df[5].astype(str)+'('+df[6]+')'
    tlocus = tloc.str.replace(" ","")
    tse = df[7].astype(str)+'-'+df[8].astype(str)
    my_data = {'tncRNA_Type': df['category'], 'tRNA': tlocus, 'tncRNA_Position_on_tRNA': tse, 'Strand': df[1], 'Sequence': df[9], 'Length': df[10].astype(str), 'Read_count' :df[0]}
    my = pd.DataFrame(my_data)
    my1 = my[(my['tncRNA_Type']=="tRF-5")|(my['tncRNA_Type']=="tRF-3")|(my['tncRNA_Type']=="tRF-3-CCA")|(my['tncRNA_Type']=="other-tRF")|(my['tncRNA_Type']=="5'tRH")|(my['tncRNA_Type']=="3'tRH-CCA")|(my['tncRNA_Type']=="3'tRH")]
    my1 = my1[my1['Read_count']>=s_cut]
    if my1.empty == False:
        my1.to_csv(file4, index= False, quoting=3, sep = '\t', escapechar=' ')

else:
    pass

#tRF-1 & leader-tRFs
file5 = file2b
file5a = file1+'.4tmp'
file5b = file1+'.5tmp'

if os.path.isfile(file5) and os.path.getsize(file5)!=0:
    dt = pd.read_csv(file5, header=None, sep='\t')
    dt[4]= dt[4].str.replace("-|\\(", "\t")
    dt[4]= dt[4].str.replace("\\)", "")
    dt[4]= dt[4].str.replace("#", "-")
    dt.to_csv(file5a, header= False, index= False, quoting=3, sep = '\t', escapechar=' ')
    
    df = pd.read_csv(file5a, header=None, sep='\t')
    nuu3 = '50'
    df.insert(11, "tRNA_len", nuu3, True)
    #tRF-1
    df.loc[(df[2].str.contains('tl_')) & (df[7]==1) & (df[10]<=50) & (df[10]>=14) & (df[0]>=s_cut), 'category'] = 'tRF-1'

    #leader-tRFs
    df.loc[(df[2].str.contains('ld_')) & (df[8]==50) & (df[10]<=50) & (df[10]>=14) & (df[0]>=s_cut), 'category'] = 'leader-tRF'
        
    df = df.loc[df['category'].str.contains('tRF-1|leader-tRF',na=False)]
    tloc = df[2]+'::'+df[3]+':'+df[4].astype(str)+'-'+df[5].astype(str)+'('+df[6]+')'
    tlocus = tloc.str.replace(" ","")
    tse = df[7].astype(str)+'-'+df[8].astype(str)
    my_data2 = {'tncRNA_Type': df['category'], 'tRNA': tlocus, 'tncRNA_Position_on_tRNA': tse, 'Strand': df[1], 'Sequence': df[9], 'Length': df[10].astype(str), 'Read_count' :df[0]}
    my_ = pd.DataFrame(my_data2)
    if my_.empty == False: #check if dataframe is not empty
        my_.to_csv(file5b, index= False, quoting=3, sep = '\t', escapechar=' ')
    else:
        pass
else:
    pass

##append both into one , and sort
f = open(sys.argv[3],'r')
nn = f.read()
pmf=int(nn)/10**6

nu_4 = file4
nuu_4 = file5b

if os.path.isfile(nu_4)==True and os.path.isfile(nuu_4)==True:
    df = pd.read_csv(nu_4, sep='\t',header=None, skiprows=1)
    df2 = pd.read_csv(nuu_4, sep='\t',header=None, skiprows=1)
    mx = df.append(df2)
    mx.sort_values([0, 1], axis=0, ascending=True, inplace=True)
    mx.columns = ["tncRNA_Type", "tRNA_info", "Position(start-end)", "Strand", "Sequence", "Length", "Read_count"]
    rpm=mx['Read_count']/pmf
    mx.insert(7, "RPM", rpm, True)
    mx.to_csv(sys.argv[4], index=False, quoting=3, sep = '\t', escapechar=' ')
elif os.path.isfile(nu_4)==True and os.path.isfile(nuu_4)==False:
    df = pd.read_csv(nu_4, sep='\t',header=None, skiprows=1)
    df.sort_values([0, 1], axis=0, ascending=True, inplace=True)
    df.columns = ["tncRNA_Type", "tRNA_info", "Position(start-end)", "Strand", "Sequence", "Length", "Read_count"]
    rpm=df['Read_count']/pmf
    df.insert(7, "RPM", rpm, True)
    df.to_csv(sys.argv[4], index=False, quoting=3, sep = '\t', escapechar=' ')
elif os.path.isfile(nu_4)==False and os.path.isfile(nuu_4)==True:
    df2 = pd.read_csv(nuu_4, sep='\t',header=None, skiprows=1)
    df2.sort_values([0, 1], axis=0, ascending=True, inplace=True)
    df2.columns = ["tncRNA_Type", "tRNA_info", "Position(start-end)", "Strand", "Sequence", "Length", "Read_count"]
    rpm=df2['Read_count']/pmf
    df2.insert(7, "RPM", rpm, True)
    df2.to_csv(sys.argv[4], index=False, quoting=3, sep = '\t', escapechar=' ')
else:
    print('No tncRNA')
