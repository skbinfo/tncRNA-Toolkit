<!DOCTYPE html>
<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
<meta name="keywords" content="tRNA, tRNA modification, tRFs, tRNA derived fragments, tncRNA, non-coding, NIPGR" />
<title>tncRNA: A pipeline for identification of tRNA-derived ncRNAs</title>
<link rel="stylesheet" href="styles.css" type="text/css" />
  <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.3.1/jquery.min.js"></script-->
  <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js"></script-->
  <script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.7.1/jquery.min.js"></script>
  <script type="text/javascript" src="jquery.dropotron-1.0.js"></script>
  <script type="text/javascript" src="jquery.slidertron-1.1.js"></script>
<style>
.topnav {
  overflow: hidden;
  background-color: #333;
  margin-top: 4%;
}

.topnav a {
  float: left;
  color: #f2f2f2;
  text-align: center;
  padding: 1% 2%;
  text-decoration: none;
  font-size: 90%;
}

.topnav a:hover {
  background-color: #759090;
}
</style>
</head>
<body>
<div class="wrapper">
<div style='display:table;width:100%;background:#7a7b76;'>
<a href="#"><img style='width:10%;background-color:white;float:left;border-radius:50%;border:solid 10px darkgrey;margin:0.2%' src=images/logo.png></a>
<div class="topnav">
				<!-- MENU -->
				<a style='color:white' href="index.php">Home</a>
                              	<a style='color:white' href="manual.php">Manual</a>
				<a style='color:white' href="data_download.php">Data download</a>
				<a style='color:white' href="team.php">Team</a>
				<a style='color:white' href="contact.php">Contact</a>
				<!-- END MENU -->
</div></div>
<br>
<button onclick="topFunction()" id="myBtn" title="Go to top">Top</button>
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN">
<html>
<head>
	<style>
td {
text-align:center;
background: lightgrey;
width:50%;
}
td > a{
color: darkblue;
font-weight:bold;
}
.center {
  display: block;
  margin-left: auto;
  margin-right: auto;
  //width: 50%;
}
	</style>
</head>
<body>
<h3>tncRNA</h3>
<!--<a style='float:right;color:blue;font-size:20px;background:yellow;border:solid 1px;' href="manual/manual.pdf" target="_blank">&#8608; PDF</a><br>-->

<table>
	<tr>
		<td><a href="#a1">1. Brief Introduction</a></td>
		<td><a href="#a2">2. Workflow of tncRNA</a></td>
	</tr>
	<tr>
		<td><a href="#a3">3. Download</a></td>
		<td><a href="#a4">4. Prerequisites</a></td>
	</tr>
	<tr>
		<td><a href="#a5">5. Installation</a></td>
		<td><a href="#a6">6. Inside the Package</a></td>
	</tr>
	<tr>
		<td><a href="#a7">7. Usage</a></td>
		 <td><a href="#a8">8. Output</a></td>
	</tr>
</table>
<br>
<div style='padding-left:5%;padding-right:5%'>
<h4 id="a1">1. Brief Introduction</h4><br>
<p align="justify">This pipeline, <b>tncRNA</b>, is designed for the identification of 
tRNA-derived small ncRNAs (tncRNAs) from high throughput sequencing data.
 This is built in the <i>python3</i>, alongwith shell/bash programming language.
 It can detect tncRNAs generating from either mature tRNAs or their leader or trailer sequences
 (upto 50 bp) viz. tRF-5, tRF-3(CCA), tRF-1, leader tRFs, 5’ tRH, 3’ tRH(CCA) 
and other-tRFs depending upon their site of cleavage on their parental tRNAs and length.</p>

<h4 id="a2">2. Worklow of tncRNA</h4><br>
 <figure>
  <img src="manual/manual_html_94adcef914d73139.png" class="center";>
  <figcaption><b>Figure 1.</b> The outline of tncRNA.</figcaption>
</figure>
<br>

<h4 id="a3">3. Download</h4>
<p align="justify">Download the tncRNA package from below link:<br>
<a href="http://172.16.2.22/tncRNA/tncRNA.centos.tar.gz" style='color:blue'>tncRNA.centos.tar.gz</a><br>
<a href="http://172.16.2.22/tncRNA/tncRNA.ubuntu.tar.gz" style='color:blue'>tncRNA.ubuntu.tar.gz</a><br>
<!-- 223.31.159.8 -->
</p>

<h4 id="a4">4. Prerequisites</h4><br>
<p align="justify">1. <a target="_blank" href="https://www.python.org/downloads/" style='color:blue'>python3</a><br>
 Python modules: pandas (v1.1.0), biopython (v1.77)<br>
 Python module can be easily installed by following command:<br>
<code>pip3 install &lt;module name&gt; --user</code><br>
2. <a target="_blank" href="http://lowelab.ucsc.edu/tRNAscan-SE/" style='color:blue'>tRNAscan-SE</a> (v2.0.6)<br>
3. <a target="_blank" href="http://www.htslib.org/download/" style='color:blue'>samtools</a> (v1.10)<br>
4. <a target="_blank" href="https://bedtools.readthedocs.io/en/latest/content/installation.html" style='color:blue'>bedtools</a> (v2.29.2)<br>
[<b>Note:</b> <i>python3</i>, <i>tRNAscan-SE</i>, <i>samtools</i>, and <i>bedtools</i> are needed to be globaly installed or included it in the path.]<br>
5. <a target="_blank" href="http://bowtie-bio.sourceforge.net/index.shtml" style='color:blue'>Bowtie1</a> (v1.3.0)<br>
6. <a target="_blank" href="http://tesla.pcbi.upenn.edu/hamr/" style='color:blue'>HAMR</a> (v1.2)<br>
[<b>Note:</b> <i>bowtie1</i>, and <i>HAMR</i> are already provided in tar package.]<br>

<h4 id="a5">5. Installation</h4><br>
Extract the tarball using the command:<br><br>
<code>tar -xf tncRNA.tar.gz</code><br><br>
<code>cd tncRNA/</code><br><br>
Installing bowtie1:<br><br>
<code>cd util/</code><br><br>
<code>unzip bowtie-1.3.0-src.zip</code><br><br>
<code>cd bowtie-1.3.0-src/</code><br><br>
<code>make</code><br><br>
<code>cd ../../</code><br><br>
[<b>Note:</b> HAMR-1.2 is already extraced in 'util/HAMR-1.2' and it required <i>python2</i> for running. If you have <i>python2</i>, no need of further installation for HAMR.<br><br>

<h4 id="a6">6. Inside the Package</h4><br>
<p align="justify">This distribution includes the python3 script 
<i>tncRNAs.py</i> and other scripts in the ‘util’ folder which 
are needed to run the main script.</p>
<figure>
  <img src="manual/ls.png" style="width:20%;">
  <figcaption><b>Figure 1.</b> Listing tncRNA directory.</figcaption>
</figure>
<br>

<h4 id="a7">7. Usage</h4><br>
First, create bowtie index for genome:<br><br>
<code>python3 tncRNAs.py -g &lt;genome fasta&gt; -s &lt;species name&gt;</code><br><br>
It will automatically create the bowtie index alongwith needed files in "lib/indexes/&lt;provided species name&gt;"
<br><br>
<b>Note</b>:<br> Genome fasta header should start with '&gt;chr[Num]'. Mitochondrial and plastid fasta headers also should be as chrMt &amp; chrPt respectively. This will be helpful for automation of scripts, and separation of nuclear &amp; organellar region.<br><br>
Once dealt with index build, user can further analyse the processed sRNA single-end data for that species.<br><br>

tncRNAs prediction:<br><br>
<p align="justify"><code>python3 tncRNAs.py -s &lt;species name&gt; -i &lt;processed small RNA reads&gt; -o &lt;output dir&gt;</code><br>
<b>Options-</b><br>
<b style='padding-left:3%'>-h</b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;print help<br>
<b style='padding-left:3%'>-s</b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;species name<br>
<b style='padding-left:3%'>-i</b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;small RNA reads (quality filtered and adapter trimmed)<br>
<font style="padding-left:6%">fastq (.fastq/.fq) file</font><br>

<b style='padding-left:3%'>-o</b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;output directory<br>

<b>Miscellaneous options-</b><br>
<b style='padding-left:3%'>-v</b>&nbsp;&nbsp;&nbsp;&lt;int&gt;&nbsp;&nbsp;&nbsp;mismatch or gap allowed [default: 0]<br>
<font style="padding-left:9%">This option is for providing the integer value for -v option to the Bowtie aligner. It is recommended not to provide value &gt;3.</font><br>
<b style='padding-left:3%'>-m</b>&nbsp;&nbsp;&nbsp;&lt;int&gt;&nbsp;&nbsp;limit to suppress all alignments if more than <int> exist [default: 50]<br>
<b style='padding-left:3%'>-c</b>&nbsp;&nbsp;&nbsp;&lt;int&gt;&nbsp;&nbsp;&nbsp;cut-off value for read count [default: 10]<br>
<b style='padding-left:3%'>-t</b>&nbsp;&nbsp;&nbsp;&lt;int&gt;&nbsp;&nbsp;&nbsp;number of threads [default: 1]<br>
</p>


<h4 id="a8">8. tncRNA output</h4><br>
<p align="justify">In the provided directory with ‘-o’ option, <b>tncRNA</b> provides three result files:<br>
&nbsp;&nbsp;&nbsp;1. <u><i>tncRNAs.csv</i></u><br/>
<font style="padding-left:3%">tncRNAs findings and their respective parental tRNA information, strand, sequence, length, read count, RPM value and Modification are stored in <i>tncRNAs.csv</i>.</font><br>
<font style="padding-left:3%">Total of nine columns in <i>tncRNAs.csv</i> (shown in Figure 3) are as follows:</font><br>

<font style="padding-left:5%">1. <i>tncRNA_Type:</i></font><br>
<font style="padding-left:6%">Different type of tncRNAs: tRF-5', tRF-3'(CCA), 3’ tRH(CCA), 5’ tRH, tRF-1, leader-tRF and other-tRF.</font><br>

<font style="padding-left:5%">2. <i>tRNA_info:</i></font><br>
<font style="padding-left:6%">This column contains tRNA information in the following format:</font><br>
<font style="padding-left:7%" color="darkblue">tRNA:amino_acid:anticodon::chromosome:start-end(strand)</font><br>
<font style="padding-left:6%">'ld_tRNA' means the leader sequence of the corresponding tRNA.</font><br>
<font style="padding-left:6%">It is 50 bp upstream to the chromosome locus of respective tRNA.</font><br>
<font style="padding-left:6%">While, 'tl_tRNA' is the trailer sequence, 50 bp downstream to tRNA sequence.</font><br>

<font style="padding-left:5%">3. <i>Position(start-end):</i></font><br>
<font style="padding-left:6%">This is the aligned-position of tncRNA on its respective origin 
(mature tRNA, trailer or leader sequence).</font><br>

<font style="padding-left:5%">4. <i>Strand:</i></font><br>
<font style="padding-left:6%">The +ve/-ve strand of tRNA origin, on which tncRNAs are aligned.</font><br>

<font style="padding-left:5%">5. <i>Sequence:</i> tncRNA sequence </font><br>

<font style="padding-left:5%">6. <i>Length:</i> Length of tncRNA</font><br>

<font style="padding-left:5%">7. <i>Read_count:</i></font><br>
<font style="padding-left:6%">Total number of mapped reads from small RNA-Seq supporting the tncRNA.</font><br>

<font style="padding-left:5%">8. <i>RPM:</i> Reads per million</font><br>
<font style="padding-left:5%">9. <i>Modification:</i> [base position on tRNA]base-&gt;Modification</font><br>
<font style="padding-left:6%">"NF" in column, if no modification found for tncRNAs.</font><br>
</p>
<figure><img src="manual/manual_html_3a0cedf558dfe5be.png" style="width:90%;"><figcaption><b>Figure 3.</b>
 Output file <i>tncRNAs.csv</i> generated from <b>tncRNA</b> showing relevant information about each identified tncRNA.</figcaption></figure>
<br/>

2. <u>Alignment.txt</u><br>
<p align="justify">This file depicts the alignment of each tncRNA sequence over
 its origin sequences, <i>i.e.</i>, mature tRNA, 5’ leader and 3’ trailer precursor 
tRNA sequences, including basic meta-information from <i>tncRNAs.csv</i> (Figure 4). 
Scoring is provided as per parameters:<br/>
<i>identical=1, non-Identical=-1, gap-open=-1, gap-extend=-0.5</i></p>
<figure><img src="manual/manual_html_dd56ff91c80b6639.png" style="width:50%;"><figcaption>
<b>Figure 4.</b> Alignment file <i>Alignment.txt</i> showing the alignment 
of tncRNA on its respective parental tRNA.</figcaption></figure>
<br/>

3. <u>Result.log</u><br/>
<p align="justify">The statistical information regarding the total count of different tncRNA 
sub-types are stored by this log file (Figure 5).</p>
<figure><img src="manual/manual_html_795c376f1bac2b77.png" style="width:12%;"><figcaption>
<b>Figure 5.</b> The log file <i>Result.log</i> generated by <b>tncRNA</b> showing 
statistical information for each sample analyzed.</figcaption></figure><br/>
<p align="justify"><u>Note:</u> ‘tncRNAs.csv’ and ‘Alignment.txt’ files will not 
be generated if no tncRNAs are detected in a sample!</p>
<br/>
<br>
</div>
</body>
</html>
<script>
// When the user scrolls down 20px from the top of the document, show the button
window.onscroll = function() {scrollFunction()};

function scrollFunction() {
  if (document.body.scrollTop > 20 || document.documentElement.scrollTop > 20) {
    document.getElementById("myBtn").style.display = "block";
  } else {
    document.getElementById("myBtn").style.display = "none";
  }
}

// When the user clicks on the button, scroll to the top of the document
function topFunction() {
  document.body.scrollTop = 0;
  document.documentElement.scrollTop = 0;
}
</script>

<div class="clear">			
<p class="footer">
<a href="https://scholar.google.co.in/citations?user=6sFjjocAAAAJ&hl=en" target="_blank">Dr. Shailesh Kumar</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
<a href="http://www.nipgr.res.in/home/home.php" target="_blank">NIPGR</a> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
<a href="http://www.dbtindia.nic.in/" target="_blank">DBT</a> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 
<a href="http://www.nipgr.res.in/research/dr_shailesh.php" target="_blank">&copy; SKLab</a>
</p>
</div>
</body>
</html>
