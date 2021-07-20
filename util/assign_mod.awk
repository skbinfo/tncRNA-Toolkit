#!/usr/bin/awk -f
#Assign the modification to tncRNA
BEGIN{FS="\t";}
{
	c=0
	x=split($3,array1,"-");
	z=split($9,array4,";|]");
}
{if(z>1){
	printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t",$1,$2,$3,$4,$5,$6,$7,$8;
	for( i = 1; i <= z; i+=2 )
		if(array4[i]>=array1[1] && array4[i]<=array1[2])
		printf "%s;","["array4[i]"]"array4[c+=2];
	printf "\n";
}
else if(z==1){print $0}
}
OFS="\t"
