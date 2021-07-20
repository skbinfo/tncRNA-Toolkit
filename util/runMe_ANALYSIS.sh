#! /bin/bash

#command args 
while getopts 'x:i:v:m:o:l:t:p:h' c
do
	case $c in
		x) index="$OPTARG";;
		i) reads="$OPTARG";;
		v) mismatches="$OPTARG";;
		m) multimap="$OPTARG";;
		o) outdir="$OPTARG";;
		t) threads="$OPTARG";;
		l) label="$OPTARG";;
		p) path="$OPTARG";;
		h) echo ""
		   echo "   -x  <path to dir created by 'runMe_INDEX.sh'>"
		   echo "   -i  <fastq file>"
		   echo	"   -o  <outdir file name>"
		   echo	"   -l  <set label for intermediate file>"
		   echo ""
		   echo "  Miscellenious:"
		   echo "   -t  <int> number of threads [default:1]"
		   echo "   -v  <int> report end-to-end hits w/ <=v mismatches [default: 0]"
		   echo "   -m  <int> suppress all alignments if > <int> exist [default: 50]"
		   echo ""
		   exit 1;;
	   [?]) printf "\n      usage: runMe_ANALYSIS.sh -x [index dir] -i [fastq reads] -o [output dir] -l [label for files]\n\n"
		   exit 1;;	
	esac
done
shift $((OPTIND-1))

#check command & tools

bpath="${path}/util/bowtie-1.3.0-src"
hamr_path="${path}/util/HAMR-1.2"
hamr="python2 ${hamr_path}/hamr.py"
hamr_model="${hamr_path}/models/euk_trna_mods.Rdata"


if ! command -v python2 &> /dev/null
then
        echo "python2 not found"
        exit 1;
elif ! command -v ${hamr_path}/hamr.py &> /dev/null
then
        echo "${hamr_path}/hamr.py not found or not in path"
        exit 1;
fi

if [[ $reads =~ .fastq$ ]];then
	bn=`basename $reads .fastq`
elif [[ $reads =~ .fq$ ]];then
	bn=`basename $reads .fq`
fi

#SAM: column 2: the reads with flag 4 are unmapped, the reads with flag 0 are mapped to the forward strand and the reads with flag 16 are mapped to the reverse strand.

echo "[`date +%d-%m-%Y\ %H:%M:%S`] Aligning the filtered reads to bowtie index ... "
if [ -z $threads ];then
	${bpath}/bowtie -v 2 --best -x ${index}/index/bt_index $reads -S ${outdir}/${bn}.sam
elif [ ! -z $threads ];then
	${bpath}/bowtie -v 2 --best -p $threads -x ${index}/index/bt_index $reads -S ${outdir}/${bn}.sam
fi

delv="${outdir}/${bn}.sam"

echo "[`date +%d-%m-%Y\ %H:%M:%S`] Analyse for modification sites... "
samtools sort ${outdir}/${bn}.sam -o ${outdir}/${bn}.sort.sam
samtools view -bS ${outdir}/${bn}.sort.sam -o ${outdir}/${bn}.bam
$hamr ${outdir}/${bn}.bam ${index}/genome_artificial.fa ${hamr_model} ${outdir}/mod_res ${bn} 30 10 0.05 H4 0.01 0.05 0.05 1>>${outdir}/hamr.run 2>>${outdir}/hamr.run

delv+=" ${outdir}/${bn}.sort.sam ${outdir}/${bn}.bam ${outdir}/mod_res"

MODDIR=${outdir}/mod_res
MODFILE=${outdir}/mod_res/${bn}.mods.txt
if [ -d "$MODDIR" ]; then
	if [ -f "$MODFILE" ]; then
		grep "^tRNA\|^ld_\|^tl_" ${outdir}/mod_res/${bn}.mods.txt|grep -wP "\bTRUE\b" \
			|awk '{print $1"\t"$2+1"]"$4"->"$NF}' \
			|awk '{if(a[$1])a[$1]=a[$1]";"$2; else a[$1]=$2;}END{for (i in a)print i, a[i];}' OFS='\t' \
			|sed '1 i tRNA_info\tModification' > ${outdir}/tRNA_modification.txt
		rm -f ${outdir}/hamr.run
		find ${outdir}/tRNA_modification.txt -empty -exec rm -f "{}" \;
		rec=${outdir}/tRNA_modification.txt
		if [ ! -f "$rec" ]; then
			echo "-- No Modification site found on mature tRNAs/leader-tRNAs/trailer-tRNAs region --"
		fi
	else
		echo ""
		echo "-- No Modification site found on mature tRNAs/leader-tRNAs/trailer-tRNAs region --"
		echo ""
	fi
else
	echo ""
	echo "-- ERROR!!! HAMR is not running --"
fi

echo "[`date +%d-%m-%Y\ %H:%M:%S`] Counting of mapped reads and unique read count from SAM ... "
awk -F\\t '!/@/{if($2==0||$2==16)print $10}' ${outdir}/${bn}.sort.sam|sort|uniq -c|sort -n|awk '{print $2"\t"$1}' \
	>${outdir}/${bn}.unique_reads.txt
TMR=`awk '{sum += $2} END {print sum}' ${outdir}/${bn}.unique_reads.txt|tr '\n' ' '|sed 's/ //g'`
printf "$TMR" > ${outdir}/Total_mapped_reads
echo "[`date +%d-%m-%Y\ %H:%M:%S`] Converting of unique read to fasta format ... "
awk '{print ">"$2"\n"$1}' ${outdir}/${bn}.unique_reads.txt > ${outdir}/${bn}.unique_reads.fa
delv+=" ${outdir}/${bn}.unique_reads.txt ${outdir}/${bn}.unique_reads.fa"

echo "[`date +%d-%m-%Y\ %H:%M:%S`] Aligning of mapped unique reads ... "
if [ -z $threads ];then
	if [ -z $mismatches ] && [ -z $multimap ];then
		${bpath}/bowtie -f --norc -v 0 -m 50 --best --strata -x ${index}/index/bt_index ${outdir}/${bn}.unique_reads.fa \
			| grep 'tRNA\|ld_\|tl_' > ${outdir}/bt.tab.out
	elif [ ! -z $mismatches ] && [ -z $multimap ];then
		${bpath}/bowtie -f --norc -v ${mismatches} -m 50 --best --strata -x ${index}/index/bt_index ${outdir}/${bn}.unique_reads.fa \
			| grep 'tRNA\|ld_\|tl_' > ${outdir}/bt.tab.out
	elif [ -z $mismatches ] && [ ! -z $multimap ];then
		${bpath}/bowtie -f --norc -v 0 -m ${multimap} --best --strata -x ${index}/index/bt_index ${outdir}/${bn}.unique_reads.fa \
			| grep 'tRNA\|ld_\|tl_' > ${outdir}/bt.tab.out
	elif [ ! -z $mismatches ] && [ ! -z $multimap ];then
		${bpath}/bowtie -f --norc -v ${mismatches} -m ${multimap} --best --strata -x ${index}/index/bt_index ${outdir}/${bn}.unique_reads.fa \
			| grep 'tRNA\|ld_\|tl_' > ${outdir}/bt.tab.out
	fi
else
	if [ -z $mismatches ] && [ -z $multimap ];then
		${bpath}/bowtie -f --norc -v 0 -m 50 --best --strata -p $threads -x ${index}/index/bt_index ${outdir}/${bn}.unique_reads.fa \
			| grep 'tRNA\|ld_\|tl_' > ${outdir}/bt.tab.out
	elif [ ! -z $mismatches ] && [ -z $multimap ];then
		${bpath}/bowtie -f --norc -v ${mismatches} -m 50 --best --strata -p $threads -x ${index}/index/bt_index ${outdir}/${bn}.unique_reads.fa \
                        | grep 'tRNA\|ld_\|tl_' > ${outdir}/bt.tab.out
	elif [ -z $mismatches ] && [ ! -z $multimap ];then
		${bpath}/bowtie -f --norc -v 0 -m ${multimap} --best --strata -p $threads -x ${index}/index/bt_index ${outdir}/${bn}.unique_reads.fa \
			| grep 'tRNA\|ld_\|tl_' > ${outdir}/bt.tab.out
	elif [ ! -z $mismatches ] && [ ! -z $multimap ];then
		${bpath}/bowtie -f --norc -v ${mismatches} -m ${multimap} --best --strata -p $threads -x ${index}/index/bt_index ${outdir}/${bn}.unique_reads.fa \
			| grep 'tRNA\|ld_\|tl_' > ${outdir}/bt.tab.out
	fi
fi

rm -rf $delv
