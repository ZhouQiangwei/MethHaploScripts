## This is a script to calculate the haplotype results in different way.
##

cd WGBS/imr90/
methsnpp="snpwithmr.p.hapcut2.haplo.txt"
methsnpn="snpwithmr.n.hapcut2.haplo.txt"
methp="mr.hapcut2.p.haplo.txt"
methn="mr.hapcut2.n.haplo.txt"
snphap="snp.hapcut2.haplo.txt"
## meth + snp
methsnpmerge="snpwithmr.haplo.merge.txt"
methmerge="mr.haplo.merge.txt"
hichap="/public/home/qwzhou/project/methhap/hic/imr90/hic_out/bowtie_results/bwt2/haplo.txt"

awk -v OFS="\t" '$1!~/\*/{if($1~/BLOCK/){if(chr!="" && end>start && start!=0){print "=== "chr,start,end"\n"inf};chr="";start=0;end=0;inf="";}else{if($2=="-"){if(chr!="" && end>start && start!=0){print "=== "chr,start,end"\n"inf};chr="";start=0;end=0;inf=""}else if($2==0 || $2==1){if(start==0){chr=$4;start=$5;inf=$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}else{inf=inf"\n"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7;if($5>end){end=$5}}}}}' ${methsnpp} ${methsnpp} | awk -v OFS="\t" '$1~/^=/{print $2,$3,$4}' | sort -k1,1 -k2,2n -k3,3n | bedtools merge > ${methsnpmerge}

awk '{a+=($3-$2);b++}END{print a,b,a/b}' ${methsnpmerge}
awk '{if($3-$2>2000){a++}else if($3-$2>1000){b++}else if($3-$2>500){c++}}END{print a,b,c}' ${methsnpmerge}

## meth
awk -v OFS="\t" '$1!~/\*/{if($1~/BLOCK/){if(chr!="" && end>start && start!=0){print "=== "chr,start,end"\n"inf};chr="";start=0;end=0;inf="";}else{if($2=="-"){if(chr!="" && end>start && start!=0){print "=== "chr,start,end"\n"inf};chr="";start=0;end=0;inf=""}else if($2==0 || $2==1){if(start==0){chr=$4;start=$5;inf=$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}else{inf=inf"\n"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7;if($5>end){end=$5}}}}}' ${methp} ${methp}  | awk -v OFS="\t" '$1~/^=/{print $2,$3,$4}' | sort -k1,1 -k2,2n -k3,3n | bedtools merge > ${methmerge}
awk '{if($3-$2>2000){a++}else if($3-$2>1000){b++}else if($3-$2>500){c++}}END{print a,b,c}' ${methmerge}
awk '{a+=($3-$2);b++}END{print a,b,a/b}' ${methmerge}

## snp
awk -v OFS="\t" '$1!~/\*/{if($1~/BLOCK/){if(chr!="" && end>start && start!=0){print "=== "chr,start,end"\n"inf};chr="";start=0;end=0;inf="";}else{if($2=="-"){if(chr!="" && end>start && start!=0){print "=== "chr,start,end"\n"inf};chr="";start=0;end=0;inf=""}else if($2==0 || $2==1){if(start==0){chr=$4;start=$5;inf=$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}else{inf=inf"\n"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7;if($5>end){end=$5}}}}}' ${snphap} | awk -v OFS="\t" '$1~/^=/{print $2,$3,$4}' | awk '{a+=($3-$2);b++}END{print a,b,a/b}'

awk -v OFS="\t" '$1!~/\*/{if($1~/BLOCK/){if(chr!="" && end>start && start!=0){print "=== "chr,start,end"\n"inf};chr="";start=0;end=0;inf="";}else{if($2=="-"){if(chr!="" && end>start && start!=0){print "=== "chr,start,end"\n"inf};chr="";start=0;end=0;inf=""}else if($2==0 || $2==1){if(start==0){chr=$4;start=$5;inf=$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}else{inf=inf"\n"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7;if($5>end){end=$5}}}}}' ${snphap} | awk -v OFS="\t" '$1~/^=/{print $2,$3,$4}' | awk '{if($3-$2>2000){a++}else if($3-$2>1000){b++}else if($3-$2>500){c++}}END{print a,b,c}'

## snp meth +hic
##merge snpmeth
~/project/methhap/WGBS/bsmerge ${methsnpp} ${methsnpn} > snpwithmr.merge.hapcut2.haplo.txt
##merge hic
~/project/methhap/WGBS/bsmergehic snpwithmr.merge.hapcut2.haplo.txt ${hichap} > snpwithmr.hap.mergedhic.txt

awk -v OFS="\t" 'BEGIN{a=0;b=0;chr=""}{
     if($1~/^B/){
         if(b>a){print chr,a,b;a=0;b=0;chr="";}
         else{a=0;b=0;chr="";}
     }
     else if(a==0){
         chr=$1;a=$2;b=$2
     }else{
         if($2>b+3000){
             if(b>a){print chr,a,b; a=$2;b=$2;chr=$1}
         }else if($2>b){
             b=$2
         }
     }
}END{if(b>a){print chr,a,b}}' snpwithmr.hap.mergedhic.txt | awk '{a+=($3-$2);b++}END{print a,b,a/b}'

awk -v OFS="\t" 'BEGIN{a=0;b=0;chr=""}{
     if($1~/^B/){
         if(b>a){print chr,a,b;a=0;b=0;chr="";}
         else{a=0;b=0;chr="";}
     }
     else if(a==0){
         chr=$1;a=$2;b=$2
     }else{
         if($2>b+3000){
             if(b>a){print chr,a,b; a=$2;b=$2;chr=$1}
         }else if($2>b){
             b=$2
         }
     }
}END{if(b>a){print chr,a,b}}' snpwithmr.hap.mergedhic.txt | awk '{if($3-$2>2000){a++}else if($3-$2>1000){b++}else if($3-$2>500){c++}}END{print a,b,c}'

### asm on repeat

## meth snp
bedtools intersect -a ${methsnpmerge} -b ~/project/Genome/hg38/hg38.fa.out.repeat -wo | awk '{a[$4" "$5" "$6]}END{print length(a)}'

## meth
bedtools intersect -a ${methmerge} -b ~/project/Genome/hg38/hg38.fa.out.repeat -wo | awk '{a[$4" "$5" "$6]}END{print length(a)}'

## snp
awk -v OFS="\t" '$1!~/\*/{if($1~/BLOCK/){if(chr!="" && end>start && start!=0){print "=== "chr,start,end"\n"inf};chr="";start=0;end=0;inf="";}else{if($2=="-"){if(chr!="" && end>start && start!=0){print "=== "chr,start,end"\n"inf};chr="";start=0;end=0;inf=""}else if($2==0 || $2==1){if(start==0){chr=$4;start=$5;inf=$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}else{inf=inf"\n"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7;if($5>end){end=$5}}}}}' ${snphap} | awk -v OFS="\t" '$1~/^=/{print $2,$3,$4}' | bedtools intersect -a - -b ~/project/Genome/hg38/hg38.fa.out.repeat -wo | awk '{a[$4" "$5" "$6]}END{print length(a)}'

### + hic
awk -v OFS="\t" 'BEGIN{a=0;b=0;chr=""}{
     if($1~/^B/){
         if(b>a){print chr,a,b;a=0;b=0;chr="";}
         else{a=0;b=0;chr="";}
     }
     else if(a==0){
         chr=$1;a=$2;b=$2
     }else{
         if($2>b+3000){
             if(b>a){print chr,a,b; a=$2;b=$2;chr=$1}
         }else if($2>b){
             b=$2
         }
     }
}END{if(b>a){print chr,a,b}}' snpwithmr.hap.mergedhic.txt | bedtools intersect -a - -b ~/project/Genome/hg38/hg38.fa.out.repeat -wo | awk '{a[$4" "$5" "$6]}END{print length(a)}'