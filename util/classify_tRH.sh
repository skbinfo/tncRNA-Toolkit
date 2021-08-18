#!/bin/sh
# @AJeet Singh: singh.ajeet@nipgr.ac.in
#This script reclassify the classified tncRNA based on anticodon loop
tncRNA=$1
ac_pos=$2
outdir=$3

awk '/other-tRF|tRF-1|leader-tRF/{print}' OFS='\t' $tncRNA > ${outdir}/another

head -1 $tncRNA > ${outdir}/header
sort $ac_pos > ${outdir}/ap_sort

##3 tRNA half CCA / 3 tRNA half
grep -P "3'tRH-CCA\t" $tncRNA > ${outdir}/3tRHCCA
sort -t $'\t' -k2 ${outdir}/3tRHCCA > ${outdir}/3tRHCCA_sort

join -t $'\t' -1 2 -2 1 -a1 ${outdir}/3tRHCCA_sort ${outdir}/ap_sort \
	|awk '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS='\t' \
	| awk '{gsub("-","\t",$3)}1' OFS='\t' \
	| awk '{if($3<$11 && $3>$10) print}' OFS='\t' > ${outdir}/3tRHCCA.csv

join -t $'\t' -1 2 -2 1 -a1 ${outdir}/3tRHCCA_sort ${outdir}/ap_sort \
        |awk '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS='\t' \
        | awk '{gsub("-","\t",$3)}1' OFS='\t' \
        | awk '{if($3>=$11 && $3>$10) print}' OFS='\t' \
	| awk '{print "tRF-3-CCA",$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS='\t' > ${outdir}/tRF-3CCA.csv

grep -P "3'tRH\t" $tncRNA > ${outdir}/3tRH
sort -t $'\t' -k2 ${outdir}/3tRH > ${outdir}/3tRH_sort

join -t $'\t' -1 2 -2 1 -a1 ${outdir}/3tRH_sort ${outdir}/ap_sort \
        |awk '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS='\t' \
        | awk '{gsub("-","\t",$3)}1' OFS='\t' \
        | awk '{if($3<$11 && $3>$10) print}' OFS='\t' > ${outdir}/3tRH.csv

join -t $'\t' -1 2 -2 1 -a1 ${outdir}/3tRH_sort ${outdir}/ap_sort \
        |awk '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS='\t' \
        | awk '{gsub("-","\t",$3)}1' OFS='\t' \
        | awk '{if($3>=$11 && $3>$10) print}' OFS='\t' \
        | awk '{print "tRF-3",$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS='\t' > ${outdir}/tRF-3.csv

##tRF 3 CCA / tRF 3
grep -w "tRF-3-CCA" $tncRNA > ${outdir}/tRF-3CCA
sort -t $'\t' -k2 ${outdir}/tRF-3CCA > ${outdir}/tRF-3CCA_sort

join -t $'\t' -1 2 -2 1 -a1 ${outdir}/tRF-3CCA_sort ${outdir}/ap_sort \
        |awk '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS='\t' \
        | awk '{gsub("-","\t",$3)}1' OFS='\t' \
        | awk '{if($3<$11 && $3>$10) print}' OFS='\t' >> ${outdir}/3tRHCCA.csv

join -t $'\t' -1 2 -2 1 -a1 ${outdir}/tRF-3CCA_sort ${outdir}/ap_sort \
        |awk '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS='\t' \
        | awk '{gsub("-","\t",$3)}1' OFS='\t' \
        | awk '{if($3>=$11 && $3>$10) print}' OFS='\t' \
        | awk '{print "tRF-3-CCA",$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS='\t' >> ${outdir}/tRF-3CCA.csv

grep -P "tRF-3\t" $tncRNA > ${outdir}/tRF-3
sort -t $'\t' -k2 ${outdir}/tRF-3 > ${outdir}/tRF-3_sort

join -t $'\t' -1 2 -2 1 -a1 ${outdir}/tRF-3_sort ${outdir}/ap_sort \
        |awk '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS='\t' \
        | awk '{gsub("-","\t",$3)}1' OFS='\t' \
        | awk '{if($3<$11 && $3>$10) print}' OFS='\t' >> ${outdir}/3tRH.csv

join -t $'\t' -1 2 -2 1 -a1 ${outdir}/tRF-3_sort ${outdir}/ap_sort \
        |awk '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS='\t' \
        | awk '{gsub("-","\t",$3)}1' OFS='\t' \
        | awk '{if($3>=$11 && $3>$10) print}' OFS='\t' \
        | awk '{print "tRF-3",$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS='\t' >> ${outdir}/tRF-3.csv


##5 tRNA half
grep -P "5'tRH\t" $tncRNA > ${outdir}/5tRH
sort -t $'\t' -k2 ${outdir}/5tRH > ${outdir}/5tRH_sort

join -t $'\t' -1 2 -2 1 -a1 ${outdir}/5tRH_sort ${outdir}/ap_sort \
        |awk '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS='\t' \
        | awk '{gsub("-","\t",$3)}1' OFS='\t' \
        | awk '{if($4<$11 && $4>$10) print}' OFS='\t' > ${outdir}/5tRH.csv

join -t $'\t' -1 2 -2 1 -a1 ${outdir}/5tRH_sort ${outdir}/ap_sort \
        |awk '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS='\t' \
        | awk '{gsub("-","\t",$3)}1' OFS='\t' \
        | awk '{if($4<=$10) print}' OFS='\t' \
        | awk '{print "tRF-5",$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS='\t' > ${outdir}/tRF-5.csv

##tRF 5
grep -P "tRF-5\t" $tncRNA > ${outdir}/tRF-5
sort -t $'\t' -k2 ${outdir}/tRF-5 > ${outdir}/tRF-5_sort

join -t $'\t' -1 2 -2 1 -a1 ${outdir}/tRF-5_sort ${outdir}/ap_sort \
        |awk '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS='\t' \
        | awk '{gsub("-","\t",$3)}1' OFS='\t' \
        | awk '{if($4<$11 && $4>$10) print}' OFS='\t' >> ${outdir}/5tRH.csv

join -t $'\t' -1 2 -2  1 -a1 ${outdir}/tRF-5_sort ${outdir}/ap_sort \
        |awk '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS='\t' \
        | awk '{gsub("-","\t",$3)}1' OFS='\t' \
        | awk '{if($4<=$10) print}' OFS='\t' \
        | awk '{print "tRF-5",$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS='\t' >> ${outdir}/tRF-5.csv

cat ${outdir}/3tRHCCA.csv ${outdir}/3tRH.csv ${outdir}/5tRH.csv ${outdir}/tRF-3CCA.csv ${outdir}/tRF-3.csv ${outdir}/tRF-5.csv \
       |awk '{print $1,$2,$3"-"$4,$5,$6,$7,$8,$9}' OFS='\t' | cat ${outdir}/another - |sort|cat ${outdir}/header - > ${outdir}/result

rm -f ${outdir}/3tRHCCA ${outdir}/3tRHCCA_sort ${outdir}/ap_sort \
	${outdir}/3tRH ${outdir}/3tRH_sort ${outdir}/tRF-3CCA ${outdir}/tRF-3CCA_sort \
	${outdir}/tRF-3 ${outdir}/tRF-3_sort ${outdir}/5tRH ${outdir}/5tRH_sort \
	${outdir}/tRF-5 ${outdir}/tRF-5_sort ${outdir}/another \
	${outdir}/3tRHCCA.csv ${outdir}/3tRH.csv ${outdir}/5tRH.csv ${outdir}/tRF-3CCA.csv \
       	${outdir}/tRF-3.csv ${outdir}/tRF-5.csv ${outdir}/tncRNAs.csv
mv ${outdir}/result ${outdir}/tncRNAs.csv
