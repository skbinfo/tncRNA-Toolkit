#! /usr/bin/env python3
'''
t-ncRNA: A pipeline for identification of tRNA-derived ncRNAs
v1.0 08/10/2020
Author: Ajeet Singh
email: singh.ajeet@nipgr.ac.in
'''
import sys, getopt, os, subprocess
import pandas as pd, pathlib as p
from glob import glob
from datetime import datetime
from distutils import spawn
######################################
class bcolors:
    BACC = '\033[1;33;44m'
    ENDC = '\033[0m'
class er:
    SC = '\033[0;31;47m' 
    EC = '\033[0m'
######################################
print ('\n'+bcolors.BACC+'{:-^100}'.format(' A pipeline for identification of tRNA-derived ncRNAs ')+bcolors.ENDC) #print center align
###########Funcs###########
def printhelp():
    print('\n\t\tDESCRIPTION: tncRNAs.py identifies tRNA-derived ncRNAs from single-end small RNA sequencing dataset\n')
    print('\t\tUSAGE: python3 tncRNAs.py -s <species_name> -i <read_file> -o <output_dir>\n')
    print('\t\tARGUMENTS:')
    print('\t\t-h     print help')
    print('\t\t-s     species name')
    print('\t\t-i     small RNA reads file (.fastq/.fq)')
    print('\t\t-o     output directory')
    print('\n\t\tOPTIONAL')
    print('\t\t-c <int>  cutoff value for read count [default: 10]')
    print('\t\t-t <int>  number of threads [default: 1]')
    print('\t\t-v <int>  report end-to-end hits w/ <=v mismatches [default: 0]')
    print('\t\t-m <int>  suppress all alignments if > <int> exist [default: 50]\n')   

###check excecutables####
PYTHON=spawn.find_executable("python3")
if PYTHON is None:
    print("***ERROR: python3 is not found")
    sys.exit("Please install python3 or make sure it is in the PATH")

SAMTOOLS=spawn.find_executable("samtools")
if SAMTOOLS is None:
    print("***ERROR: samtools is not found")
    sys.exit("Please install samtools 1.1x  or make sure it is in the PATH")

BEDTOOLS=spawn.find_executable("bedtools")
if BEDTOOLS is None:
    print("***ERROR: bedtools is not found")
    sys.exit("Please install bedtools 2.2x or make sure it is in the PATH")

tRNAscan=spawn.find_executable("tRNAscan-SE")
if tRNAscan is None:
    print("***ERROR: tRNAscan-SE is not found")
    sys.exit("Please install tRNAscan-SE 2x or make sure it is in the PATH")

############Define variables##############

sps_flag=False
g_flag=False
in_flag=False
out_flag=False
t_flag=False
m_flag=False
v_flag=False
MMPATH = sys.argv[0].replace('tncRNAs.py','')
if MMPATH == "":
    MMPATH = os.getcwd()
else:
    pass
MPATH=MMPATH.replace("//","/")
print(MPATH)
##set default values
sps_name=""
in_file=""
out_dir=""
m_c=str(10)
mismatches=""
multimap=""
num_threads=""
#####################################
try:
    opts, args = getopt.getopt(sys.argv[1:],"hg::s:i:v::m::o:t::c::",["help","genome","species","read","mismatches","multimap","outdir","threads","cutoff"])
except getopt.GetoptError:
	sys.exit(printhelp())
for opt, arg in opts:
    if opt == '-h':
        sys.exit(printhelp())
    elif opt in ("-s", "--species"):
        sps_name=arg
        sps_flag=True
    elif opt in ("-g", "--genome"):
        genome=arg
        g_flag=True
    elif opt in ("-i", "--read"):
        in_file=os.path.abspath(arg)
        in_flag=True
    elif opt in ("-v", "--mismatches"):
        mismatches=arg
        v_flag=True
    elif opt in ("-m", "--multimap"):
        multimap=arg
        m_flag=True
    elif opt in ("-o", "--outdir"):
        out_dir=os.path.abspath(arg)
        out_flag=True
    elif opt in ("-c","--cutoff"):
        m_c=arg
    elif opt in ("-t", "--threads"):
        num_threads=arg
        t_flag=True

if len(sys.argv) <= 1:
    sys.exit(printhelp())
##	
###Check all necessary inputs
if sps_flag==False:
	sys.exit('Please provide species name')
if sps_flag==True and g_flag==False:
    if (os.path.isdir(MPATH+'/lib/indexes/'+sps_name)):
        pass
    else:
        sys.exit('Bowtie-index are not found for '+sps_name+' in lib\nRun the following command:\n python3 tRNA_Fs.py -s <species_name> -g <genome_fasta_file>\nYou can also provide cpu threads to bowtie-build by \'-t\' option')
if sps_flag==True and g_flag==True:
    ind_arg = [MPATH+'/util/runMe_INDEX.sh','-g',genome,'-o',MPATH+'/lib/indexes/'+sps_name,'-l',sps_name,'-p',MPATH]
    try:
        subprocess.run(ind_arg, check=True)
        sys.exit
    except subprocess.CalledProcessError as e:
        sys.exit(e.stderr)
    #Delete intermediate files
    #for tmp_i in glob(MPATH+'/lib/indexes/'+sps_name+'/'+sps_name+'*'):
     #   os.remove(tmp_i)

# check if read files exists
if in_flag==True:
    if not (os.path.isfile(in_file)):
        sys.exit(er.SC+'Please check input file...Error file:'+in_file+ ' doesn\'t exist'+er.EC)
    elif os.path.getsize(in_file) == 0:
        sys.exit(er.SC+'Please check input file...Error file:'+in_file+ ' is empty'+er.EC)

if out_flag==True:
    if (os.path.isdir(out_dir)):
        sys.exit(er.SC+out_dir+' is already there. Please provide another name or remove/rename the existing one'+er.EC)
##################################################################

if in_flag==True and in_file.endswith(('.fastq', '.fq')):
    # create outdir
    args_a = ['mkdir', '-p', out_dir]
    try:
        subprocess.run(args_a, check=True)
    except subprocess.CalledProcessError as e:
        sys.exit(e.stderr)
    
    args_run = [MPATH+'/util/runMe_ANALYSIS.sh','-x',MPATH+'/lib/indexes/'+sps_name,'-i',in_file,'-v',mismatches,'-m',multimap,'-o',out_dir,'-l',sps_name,'-t',num_threads,'-p',MPATH]
    try:
        subprocess.run(args_run, check=True)
    except subprocess.CalledProcessError as e:
        sys.exit(e.stderr)
        
    # set input variable
    # 1. indexes
    TCI = MPATH+'/lib/indexes/'+sps_name+'/index/bt_index'
    merge_tab = MPATH+'/lib/indexes/'+sps_name+'/db/merge.tab'

    # set output variable
    # 1. bowtie outputs(O)
    TiO = out_dir+'/bt_aligned_hits.tmp'

    # File add sum
    inf = out_dir+'/bt.tab.out'
    outf = out_dir+'/tncRNAs.csv'
    r_c = out_dir+'/Total_mapped_reads'

    if os.path.isfile(inf) and os.path.getsize(inf) != 0:
        if os.path.isfile(r_c) and os.path.getsize(r_c) != 0:
            args_three = ['python3', MPATH+'/util/classify.py', inf, m_c, r_c, outf]
            try:
                print(datetime.now().strftime("[%d-%m-%Y %H:%M:%S]")+' Classifying tncRNAs based on cleavage position and length .....')
                subprocess.run(args_three, check=True)
                os.remove(inf)
                os.remove(out_dir+'/Total_mapped_reads')
            except subprocess.CalledProcessError as e:
                print("\nError!\n")
                sys.exit(e.stderr)
    
    if os.path.isfile(outf) and os.path.getsize(outf) != 0:
        pos = MPATH+'/lib/indexes/'+sps_name+'/anticodon_pos'
        
        re_arg = [MPATH+'/util/classify_tRH.sh', outf, pos, out_dir]
        try:
            print(datetime.now().strftime("[%d-%m-%Y %H:%M:%S]")+' Classifying tRNA halves based on anticodon loop position .....')
            subprocess.run(re_arg, check=True)
            print(datetime.now().strftime("[%d-%m-%Y %H:%M:%S]")+' Classification done')
        except subprocess.CalledProcessError as e:
            sys.exit(e.stderr)

    for tmp_s in glob(out_dir+'/*tmp'):
        os.remove(tmp_s)

    if os.path.isfile(out_dir+'/tncRNAs.csv'):
        alns = ['python3', MPATH+'/util/align.py', out_dir+'/tncRNAs.csv', merge_tab]
        try:
            print(datetime.now().strftime("[%d-%m-%Y %H:%M:%S]")+' Creating local alignment file .....')
            with open(out_dir+"/Alignment.txt", "a") as aln:
                subprocess.run(alns, check=True, stdout=aln)
            print(datetime.now().strftime("[%d-%m-%Y %H:%M:%S]")+' Local alignment done')
        except subprocess.CalledProcessError as e:
            print("\nError!\n")
    
    if os.path.isfile(out_dir+'/tncRNAs.csv') and os.path.isfile(out_dir+"/tRNA_modification.txt"):
        if os.path.getsize(out_dir+'/tncRNAs.csv') != 0 and os.path.getsize(out_dir+"/tRNA_modification.txt") != 0:
            print(datetime.now().strftime("[%d-%m-%Y %H:%M:%S]")+' Adding modification site information for tncRNAs .....')
            df1a = pd.read_csv(out_dir+"/tncRNAs.csv", sep ='\t')
            df2a = pd.read_csv(out_dir+"/tRNA_modification.txt", sep ='\t')
            df3a=df1a.merge(df2a, on='tRNA_info', how='left')
            df3a['Modification'] = df3a['Modification'].fillna('NF')
            os.remove(out_dir+'/tncRNAs.csv')
            os.remove(out_dir+'/tRNA_modification.txt')
            df3a.to_csv(out_dir+'/result', index= False, header=None, quoting=3, sep = '\t', escapechar=' ')
            mya='awk -F\'\t\' \'{if($NF=="")gsub("","NF",$NF)}1\' OFS=\'\t\''
            os.system("awk -f "+MPATH+"/util/assign_mod.awk "+out_dir+"/result | sort -u -t $'\t' -k5,5 |"+mya+"|sort |cat "+out_dir+"/header - > "+out_dir+"/tncRNAs.csv")
            os.system("sed -i '/tncRNA_Type/s/$/\tModification/g' "+out_dir+"/tncRNAs.csv")
            os.remove(out_dir+'/result')
            os.remove(out_dir+'/header')
    
    if os.path.isfile(outf) and os.path.getsize(outf) != 0:
        re2_arg = ['python3', MPATH+'/util/count.py', outf]
        try:
            fa_one = open(out_dir+"/Result.log", "a")
            subprocess.run(re2_arg, check=True, stdout=fa_one)
            fa_one.close()
            print('Results are saved in directory: '+out_dir+'\n')
        except subprocess.CalledProcessError as e:
            print("\nError!\n")
            sys.exit(e.stderr)
    
    elif not (os.path.isfile(outf)) and os.path.isfile(out_dir+"/Result.log"):
        print('No tncRNA !\nlog file saved in: '+out_dir+'/Result.log')

    else:
        pass 

if in_flag==True and not in_file.endswith(('.fastq', '.fq')):
    sys.exit("Reads file have no .fq or .fastq extension ")
############################################END#part#1#######################################
