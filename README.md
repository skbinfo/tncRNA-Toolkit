       _                    ______   _   _      _
     _| |_                 |  __  | | \ | |    / \
    |_   _|  _____   ____  | |__| | |  \| |   /   \
      | |   |  _  | |  __| |  ____| | \ | |  / /_\ \
      | |_  | | | | | |__  | | \ \  | |\  | | |   | |
      |__ / |_| |_| |____| |_|  \_\ |_| \_| |_|   |_|

                                                                      
 tncRNA: tRNA derived non-coding RNAs
         *            *   *      ***
 website: http://www.nipgr.ac.in/tncRNA
                                                                      
 Bioinformatics Laboratory #202, National Institute of Plant Genome Research (NIPGR) New Delhi, INDIA
-------------------------------------------------------------------------------------------------------

 tncRNA is a workflow to detect tRNA derived non-coding RNAs, viz. tRF-5s, tRF-3s, tRF-1, and tRNA halves.
 It is written in the Python(v3) programming language, with additional shell scripts. Input consists of 
 genome fasta and single-end small RNA reads fastq file. tncRNA gives output files in tsv format.

 For more usage and output details, visit Manual page at "http://nipgr.ac.in/tncRNA/".

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
