#!/bin/bash
# @AJeet Singh: singh.ajeet@nipgr.ac.in
#This script identifiy the tRNA; filter tRNA; mask genome; create artifical genome & build bowtie indexes
while getopts 'g:o:l:t:p:' c
do
  case $c in
    g) genome="$OPTARG";;
    o) outdir="$OPTARG";;
    l) label="$OPTARG";;
    t) threads="$OPTARG";;
    p) path="$OPTARG";;
  esac
done 
shift $((OPTIND-1))
bpath="${path}/util/bowtie-1.3.0-src"

if ! command -v samtools &> /dev/null
then
        echo "samtools not found or not in path"
        exit 1;
elif ! command -v tRNAscan-SE &> /dev/null
then
	echo "tRNAscan-SE not found or not in path"
        exit 1;
elif ! command -v ${bpath}/bowtie &> /dev/null
then
        echo "bowtie not found or not in path"
        exit 1;
elif ! command -v ${bpath}/bowtie-build &> /dev/null
then
        echo "bowtie-build not found"
        exit 1;
fi

if [ -d "$outdir" ]
then
	echo "${outdir} is already present. Provide other name";
	exit 1;
else
	mkdir -p "$outdir"
fi

if [ -f "$genome" ]
then
	:
else
	echo "$genome not exist"
	exit 1;
fi

#program paths

cD=`pwd`
d1=`date`
printf "[$d1] Separating Nuclear region and organellar region from genome ...\n"
grep '>' ${genome}|cut -f1 -d ' '|sed 's/>//g'|grep -v 'chrMt\|chrPt' > ${outdir}/nuc_list_id.txt
grep '>' ${genome}|cut -f1 -d ' '|sed 's/>//g'|grep 'chrMt\|chrPt' > ${outdir}/org_list_id.txt

xargs samtools faidx ${genome} < ${outdir}/nuc_list_id.txt > ${outdir}/genome.nuc.fa
xargs samtools faidx ${genome} < ${outdir}/org_list_id.txt > ${outdir}/genome.org.fa

delv="${outdir}/nuc_list_id.txt ${outdir}/org_list_id.txt ${outdir}/genome.nuc.fa ${outdir}/genome.org.fa"

printf "[`date`] Scanning for nuclear tRNA ...\n"
tRNAscan-SE -Q -o ${outdir}/${label}.nuc.tRNA.csv -b ${outdir}/${label}.nuc.tRNA.bed12 -f ${outdir}/${label}.nuc.tRNA.ss ${outdir}/genome.nuc.fa

delv+=" ${outdir}/${label}.nuc.tRNA.csv ${outdir}/${label}.nuc.tRNA.bed12 ${outdir}/${label}.nuc.tRNA.ss"
printf "[`date`] Scanning for organellar tRNA ...\n"
tRNAscan-SE -Q -O -o ${outdir}/${label}.org.tRNA.csv -b ${outdir}/${label}.org.tRNA.bed12 -f ${outdir}/${label}.org.tRNA.ss ${outdir}/genome.org.fa

cat ${outdir}/${label}.nuc.tRNA.ss ${outdir}/${label}.org.tRNA.ss | grep -w 'Length\|Type' |tr '\n' '\t'|sed 's/\t$/\n/g'|sed 's/\tchr/\nchr/g'|tr '[:blank:]' ' '|sed 's/).*Type: / /g'|sed 's/ Anticodon: \|).*//g'|awk '{gsub("\\(","",$2)}1'|sed 's/ (.*\| at//g'|awk -F'\\.| ' '{print $2"-"$4"::"$1,$3,$5}'|awk '{gsub("-"," ",$2);gsub("-"," ",$3)}1'|awk '{if($3>$2)print $1":"$2-1"-"$3"(+)\t"$4-3"\t"$5+3; else{print $1":"$3-1"-"$2"(-)\t"$4-3"\t"$5+3}}'|sed 's/trna/tRNA/g' > ${outdir}/anticodon_pos

delv+=" ${outdir}/${label}.org.tRNA.csv ${outdir}/${label}.org.tRNA.bed12 ${outdir}/${label}.org.tRNA.ss"
##cat csv files and exclude the pseudo and score < 50 tRNA
cat ${outdir}/${label}.nuc.tRNA.csv ${outdir}/${label}.org.tRNA.csv|grep -v "Sequence\|Name\|---\|pseudo" | \
       awk '{if($9>50)print $1".tRNA"$2}' > ${outdir}/${label}.tRNA_header

delv+=" ${outdir}/${label}.tRNA_header"
##cat bed12 files and fetch filtered entries by tRNA_header
cat ${outdir}/${label}.nuc.tRNA.bed12 ${outdir}/${label}.org.tRNA.bed12 > ${outdir}/${label}.tRNA.bed12
grep -wf ${outdir}/${label}.tRNA_header ${outdir}/${label}.tRNA.bed12 \
	|awk '{gsub(".*.tRNA","tRNA",$4)}1' OFS='\t' > ${outdir}/${label}.tRNA.filter.bed12

delv+=" ${outdir}/${label}.tRNA.bed12 ${outdir}/${label}.tRNA.filter.bed12"
##mask pre-tRNA(tRNAs with 50 nt 5' & 3' flank) regions from genome
##bed12 for 50 nt 5' flank; remove malformed entries
awk -F\\t '{gsub("tRNA","ld_tRNA",$4);print $1,$2-50,$2,$4,$5,$6}' OFS='\t' ${outdir}/${label}.tRNA.filter.bed12 | \
	sed '/\t-[0-9]/d' >${outdir}/${label}.tRNA_up.bed12

delv+=" ${outdir}/${label}.tRNA_up.bed12"
##bed12 for 50 nt 3' flank; remove malformed entries
awk -F\\t '{gsub("tRNA","tl_tRNA",$4);print $1,$3,$3+50,$4,$5,$6}' OFS='\t' ${outdir}/${label}.tRNA.filter.bed12 \
	> ${outdir}/${label}.tRNA_down.tmp.bed12

delv+=" ${outdir}/${label}.tRNA_down.tmp.bed12"

cut -f1,2 ${genome}.fai >${outdir}/chr_len_genome.tmp #get fasta length from fai index
delv+=" ${outdir}/chr_len_genome.tmp"

#add chr length to each entry; if length exceed, change to malformed; remove malformed

sort -k1,1 ${outdir}/${label}.tRNA_down.tmp.bed12 > ${outdir}/${label}.tRNA_down.tmp.sorted.bed12
sort -k1,1 ${outdir}/chr_len_genome.tmp > ${outdir}/chr_len_genome

join -t $'\t' -1 1 -2 1 ${outdir}/${label}.tRNA_down.tmp.sorted.bed12 ${outdir}/chr_len_genome \
	| awk -F\\t '{if($3>$7)print "malformed";else{print}}' \
	|sed '/^malformed/d' |cut -f1,2,3,4,5,6 > ${outdir}/${label}.tRNA_down.bed12
delv+=" ${outdir}/chr_len_genome ${outdir}/${label}.tRNA_down.bed12 ${outdir}/${label}.tRNA_down.tmp.sorted.bed12"

cut -f1,2,3,4,5,6 ${outdir}/${label}.tRNA.filter.bed12 >${outdir}/tmp1.bed
delv+=" ${outdir}/tmp1.bed"
#get bedfile for mask the genome
cat ${outdir}/tmp1.bed ${outdir}/${label}.tRNA_up.bed12 ${outdir}/${label}.tRNA_down.bed12 >${outdir}/mask.bed12
delv+=" ${outdir}/mask.bed12"
## mask genome from pre-tRNA region
bedtools maskfasta -fi ${genome} -fo ${outdir}/genome.masked.fa -mc N -bed ${outdir}/mask.bed12
delv+=" ${outdir}/genome.masked.fa"

##remove introns, make fasta from filtered tRNA bed12, add CCA tail
bedtools getfasta -fi ${genome} -bed ${outdir}/${label}.tRNA.filter.bed12 -name -split -s -fo ${outdir}/${label}.mature_tRNA.fa
sed -i '/>/!s/$/CCA/g' ${outdir}/${label}.mature_tRNA.fa

delv+=" ${outdir}/${label}.mature_tRNA.fa"
#make fasta for 50 nt 5' tRNA flank
bedtools getfasta -name -s -fi $genome -bed ${outdir}/${label}.tRNA_up.bed12 -fo ${outdir}/${label}.up_tRNA.fa

delv+=" ${outdir}/${label}.up_tRNA.fa"
#make fasta for 50 nt 3' tRNA flank
bedtools getfasta -name -s -fi $genome -bed ${outdir}/${label}.tRNA_down.bed12 -fo ${outdir}/${label}.down_tRNA.fa

delv+=" ${outdir}/${label}.down_tRNA.fa"
##add pre-tRNAs flanks and mature tRNA as extra chromosoms to the genome (get the artificial genome)
cat ${outdir}/genome.masked.fa ${outdir}/${label}.up_tRNA.fa ${outdir}/${label}.mature_tRNA.fa ${outdir}/${label}.down_tRNA.fa \
	> ${outdir}/genome_artificial.fa

##make tab delimited fasta for pre-tRNAs
mkdir -p ${outdir}/db
cat ${outdir}/${label}.up_tRNA.fa ${outdir}/${label}.mature_tRNA.fa ${outdir}/${label}.down_tRNA.fa \
	| tr '\n' '\t'|sed 's/\t\{1,5\}$/\n/g'|sed 's/^>//g;s/\t>/\n/g' \
	|awk -F\\t '{print $1,toupper($2)}' OFS='\t'|sed '1 i tRNA_info\ttseq' > ${outdir}/db/merge.tab

##generate bowtie indexes for artificial genome 
mkdir -p ${outdir}/index
cd ${outdir}
if [ -z $threads ]
then
	${bpath}/bowtie-build genome_artificial.fa bt_index
else
	${bpath}/bowtie-build --threads $threads genome_artificial.fa bt_index
fi
mv bt_index* index

##delete intermediate
rm -rf $delv

cd ${cD}
d4=`date`
printf "\n[$d4] Complete. \n"
