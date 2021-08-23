       _                    ______   _   _      _
     _| |_                 |  __  | | \ | |    / \
    |_   _|  _____   ____  | |__| | |  \| |   /   \
      | |   |  _  | |  __| |  ____| | \ | |  / /_\ \
      | |_  | | | | | |__  | | \ \  | |\  | | |   | |
      |__ / |_| |_| |____| |_|  \_\ |_| \_| |_|   |_| Toolkit

                                                                      
 tncRNA: tRNA derived non-coding RNAs

 tncRNA Toolkit is a workflow to detect tRNA derived non-coding RNAs, viz. tRF-5s, tRF-3s, tRF-1, and tRNA halves.
 Input consists of genome fasta and single-end small RNA reads fastq file.

 For more usage and output details, visit Manual page at "http://nipgr.ac.in/tncRNA/"



[![DOI](https://zenodo.org/badge/387762887.svg)](https://zenodo.org/badge/latestdoi/387762887)


-------------------------------------------------------------------------------------------------------

## Dependencies:
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

## Usage

First, create bowtie index for genome:

<code>python3 tncRNAs.py -g < genome fasta > -s < species name ></code>

It will automatically create the bowtie index alongwith needed files in "lib/indexes/< species name >"

Note:
       
Genome fasta header should start with '>chr[Num]', for example >chr1, >chr2 so on. Mitochondrial and plastid fasta headers also should be as >chrMt & >chrPt respectively. This will be helpful for automation of scripts, and separation of nuclear & organellar region.

Once dealt with index build, user can further analyse the processed sRNA single-end data for that species.

tncRNAs prediction:

<code>python3 tncRNAs.py -s < species name > -i < processed small RNA reads > -o < output dir ></code>
       
| Options | Usage |
|---------|-------|
| -h |     print help | 
| -s |   species name |
| -i |   small RNA reads (quality filtered and adapter trimmed) fastq (.fastq/.fq) file |
| -o |   output directory |
| **Miscellaneous** | |
| -v   <int> |  mismatch or gap allowed [default: 0]. This option is for providing the integer value for -v option to the Bowtie aligner. It is recommended not to provide value >3. |
| -m   <int> | limit to suppress all alignments if more than exist [default: 50] |
| -c   <int> | cut-off value for read count [default: 10] |
| -t   <int> | number of threads [default: 1] |

--------------------------------------------------------------------------------------------------------------------
 If you have any questions, bug reports, or suggestions, please e-mail

   Dr. Shailesh Kumar
   
   Staff Scientist, Bioinformatics Laboratory #202
   
   National Institute of Plant Genome Research (NIPGR), New Delhi
    
   shailesh@nipgr.ac.in
