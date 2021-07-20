       _                    ______   _   _      _
     _| |_                 |  __  | | \ | |    / \
    |_   _|  _____   ____  | |__| | |  \| |   /   \
      | |   |  _  | |  __| |  ____| | \ | |  / /_\ \
      | |_  | | | | | |__  | | \ \  | |\  | | |   | |
      |__ / |_| |_| |____| |_|  \_\ |_| \_| |_|   |_|

                                                                      
 tncRNA: tRNA derived non-coding RNAs

 website: http://www.nipgr.ac.in/tncRNA
                                                                      
 Bioinformatics Laboratory #202, National Institute of Plant Genome Research (NIPGR) New Delhi, INDIA
-------------------------------------------------------------------------------------------------------

 tncRNA is a workflow to detect tRNA derived non-coding RNAs, viz. tRF-5s, tRF-3s, tRF-1, and tRNA halves.
 It is written in the Python(v3) programming language, with additional shell scripts. Input consists of 
 genome fasta and single-end small RNA reads fastq file. tncRNA gives output files in tsv format.

 For more usage and output details, visit Manual page at "http://nipgr.ac.in/tncRNA/".

-------------------------------------------------------------------------------------------------------
 1. Brief Introduction

This pipeline, tncRNA, is designed for the identification of tRNA-derived small ncRNAs (tncRNAs) from high throughput sequencing data. This is built in the python3, alongwith shell/bash programming language. It can detect tncRNAs generating from either mature tRNAs or their leader or trailer sequences (upto 50 bp) viz. tRF-5, tRF-3(CCA), tRF-1, leader tRFs, 5’ tRH, 3’ tRH(CCA) and other-tRFs depending upon their site of cleavage on their parental tRNAs and length.
2. Worklow of tncRNA

Figure 1. The outline of tncRNA.

3. Download

Download the tncRNA package from below link:
tncRNA.centos.tar.gz
tncRNA.ubuntu.tar.gz
4. Prerequisites

1. python3
Python modules: pandas (v1.1.0), biopython (v1.77)
Python module can be easily installed by following command:
pip3 install <module name> --user
2. tRNAscan-SE (v2.0.6)
3. samtools (v1.10)
4. bedtools (v2.29.2)
[Note: python3, tRNAscan-SE, samtools, and bedtools are needed to be globaly installed or included it in the path.]
5. Bowtie1 (v1.3.0)
6. HAMR (v1.2)
[Note: bowtie1, and HAMR are already provided in tar package.]
5. Installation

Extract the tarball using the command:

tar -xf tncRNA.tar.gz

cd tncRNA/

Installing bowtie1:

cd util/

unzip bowtie-1.3.0-src.zip

cd bowtie-1.3.0-src/

make

cd ../../

[Note: HAMR-1.2 is already extraced in 'util/HAMR-1.2' and it required python2 for running. If you have python2, no need of further installation for HAMR.

6. Inside the Package

This distribution includes the python3 script tncRNAs.py and other scripts in the ‘util’ folder which are needed to run the main script.
Figure 1. Listing tncRNA directory.

7. Usage

First, create bowtie index for genome:

python3 tncRNAs.py -g <genome fasta> -s <species name>

It will automatically create the bowtie index alongwith needed files in "lib/indexes/<provided species name>"

Note:
Genome fasta header should start with '>chr[Num]'. Mitochondrial and plastid fasta headers also should be as chrMt & chrPt respectively. This will be helpful for automation of scripts, and separation of nuclear & organellar region.

Once dealt with index build, user can further analyse the processed sRNA single-end data for that species.

tncRNAs prediction:

python3 tncRNAs.py -s <species name> -i <processed small RNA reads> -o <output dir>
Options-
-h     print help
-s     species name
-i     small RNA reads (quality filtered and adapter trimmed)
fastq (.fastq/.fq) file
-o     output directory
Miscellaneous options-
-v   <int>   mismatch or gap allowed [default: 0]
This option is for providing the integer value for -v option to the Bowtie aligner. It is recommended not to provide value >3.
-m   <int>  limit to suppress all alignments if more than exist [default: 50]
-c   <int>   cut-off value for read count [default: 10]
-t   <int>   number of threads [default: 1]
8. tncRNA output

In the provided directory with ‘-o’ option, tncRNA provides three result files:
   1. tncRNAs.csv
tncRNAs findings and their respective parental tRNA information, strand, sequence, length, read count, RPM value and Modification are stored in tncRNAs.csv.
Total of nine columns in tncRNAs.csv (shown in Figure 3) are as follows:
1. tncRNA_Type:
Different type of tncRNAs: tRF-5', tRF-3'(CCA), 3’ tRH(CCA), 5’ tRH, tRF-1, leader-tRF and other-tRF.
2. tRNA_info:
This column contains tRNA information in the following format:
tRNA:amino_acid:anticodon::chromosome:start-end(strand)
'ld_tRNA' means the leader sequence of the corresponding tRNA.
It is 50 bp upstream to the chromosome locus of respective tRNA.
While, 'tl_tRNA' is the trailer sequence, 50 bp downstream to tRNA sequence.
3. Position(start-end):
This is the aligned-position of tncRNA on its respective origin (mature tRNA, trailer or leader sequence).
4. Strand:
The +ve/-ve strand of tRNA origin, on which tncRNAs are aligned.
5. Sequence: tncRNA sequence
6. Length: Length of tncRNA
7. Read_count:
Total number of mapped reads from small RNA-Seq supporting the tncRNA.
8. RPM: Reads per million
9. Modification: [base position on tRNA]base->Modification
"NF" in column, if no modification found for tncRNAs.
Figure 3. Output file tncRNAs.csv generated from tncRNA showing relevant information about each identified tncRNA.

2. Alignment.txt

This file depicts the alignment of each tncRNA sequence over its origin sequences, i.e., mature tRNA, 5’ leader and 3’ trailer precursor tRNA sequences, including basic meta-information from tncRNAs.csv (Figure 4). Scoring is provided as per parameters:
identical=1, non-Identical=-1, gap-open=-1, gap-extend=-0.5
Figure 4. Alignment file Alignment.txt showing the alignment of tncRNA on its respective parental tRNA.

3. Result.log

The statistical information regarding the total count of different tncRNA sub-types are stored by this log file (Figure 5).
Figure 5. The log file Result.log generated by tncRNA showing statistical information for each sample analyzed.

Note: ‘tncRNAs.csv’ and ‘Alignment.txt’ files will not be generated if no tncRNAs are detected in a sample!

-------------------------------------------------------------------------------------------------------
 INSTALLATION:

 The pre-requisites includes:
 1. python3
    
    Python modules: pandas  (v1.1.0), biopython (v1.77)
    Python module can be easily installed by following command:

      pip3 install <module_name> --user

 2. tRNAscan-SE (v2.0.6)
 3. samtools (v1.10)
 4. bedtools (v2.29.2)
 Above mentioned are needed to be globally insatlled or included in the path.

 5. Bowtie1 (v1.3.0)
 6. HAMR (v1.2)
 These are provided in the tncRNA tarball.

 Extract tar package
    
     tar -xf tncRNA.centos.tar.gz
     cd tncRNA.centos/

 Install bowtie
     
     cd util/
     unzip bowtie-1.3.0-src.zip
     cd bowtie-1.3.0-src/
     make

     cd ../../

  Check HAMR is running (It needs python2). If not showing help option, please install python2
    
     python2 util/HAMR-1.2/hamr.py
  
------------------------------------------------------------------------------------------------------
 If you have any questions, bug reports, or suggestions, please e-mail

   Dr. Shailesh Kumar
   shailesh@nipgr.ac.in
