       _                    ______   _   _      _
     _| |_                 |  __  | | \ | |    / \
    |_   _|  _____   ____  | |__| | |  \| |   /   \
      | |   |  _  | |  __| |  ____| | \ | |  / /_\ \
      | |_  | | | | | |__  | | \ \  | |\  | | |   | |
      |__ / |_| |_| |____| |_|  \_\ |_| \_| |_|   |_| Toolkit

                                                                      
 tncRNA: tRNA derived non-coding RNAs
                                                                      
 Bioinformatics Laboratory 202, National Institute of Plant Genome Research (NIPGR) New Delhi, INDIA
-------------------------------------------------------------------------------------------------------

 tncRNA Toolkit is a workflow to detect tRNA derived non-coding RNAs, viz. tRF-5s, tRF-3s, tRF-1, and tRNA halves.
 Input consists of genome fasta and single-end small RNA reads fastq file.

 For more usage and output details, visit Manual page at "http://nipgr.ac.in/tncRNA/"

-------------------------------------------------------------------------------------------------------

## Dependecies:
1. python3

Python modules: pandas (v1.1.0), biopython (v1.77)

Python module can be easily installed by following command:

<code>pip3 install < module name > --user</code>

2. tRNAscan-SE (v2.0.6)
 
3. samtools (v1.10)

4. bedtools (v2.29.2)
       
[Note: python3, tRNAscan-SE, samtools, and bedtools are needed to be globaly installed or included it in the path.]
       
5. Bowtie1 (v1.3.0)
       
6. HAMR (v1.2)
       
[Note: bowtie1, and HAMR are already provided in tar package.]
       
## Installation:
Download
       
<code>git clone https://github.com/skbinfo/tncRNA-Toolkit.git</code>

<code>cd tncRNA-Toolkit/</code>

Installing bowtie1:

<code>cd util/</code>

<code>unzip bowtie-1.3.0-src.zip</code>

<code>cd bowtie-1.3.0-src/</code>

<code>make</code>

<code>cd ../../</code>

[Note: HAMR-1.2 is already extraced in 'util/HAMR-1.2' and it required python2 for running.]
4. Inside the Package

This distribution includes the python3 script tncRNAs.py and other scripts in the ‘util’ folder which are needed to run the main script.
Figure 1. Listing tncRNA directory.

5. Usage

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

Note: ‘tncRNAs.csv’ and ‘Alignment.txt’ files will not be generated if no tncRNAs are detected in a sample. 
------------------------------------------------------------------------------------------------------
 If you have any questions, bug reports, or suggestions, please e-mail

   Dr. Shailesh Kumar
   shailesh@nipgr.ac.in
