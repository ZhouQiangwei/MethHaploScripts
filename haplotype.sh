## haplotype
### haplotype assembly
### meth
~/software_devp/methyhaplo/bin/methyhaplo -M hap -a Y -m methratio.txt -s bs.align.sam -o mr.hap
### meth + snp
~/software_devp/methyhaplo/bin/methyhaplo -M hap -a Y -m methratio.txt -v snps.vcf -s bs.align.sam -o mrsnp.hap
### snp
~/software/HapCUT2/build/extractHAIRS --bam bs.align.bam --VCF snps.vcf --out snphap.fragment_file
~/software/HapCUT2/build/HAPCUT2 --fragments snphap.fragment_file --vcf snps.vcf --output snphap.haplo.txt
### hic haplotype assembly in hic result folder
samtools merge -O BAM -@ 6 k562.hic.bam ./merged/*bwt2pairs.bam
samtools sort -@ 6 -o k562.hic.sort.bam k562.hic.bam
input="k562.hic.sort.bam"  ## methyHiC.sort.p.bam
vcf="wgbs.replace.mdups.snp.filtered.sort.vcf"  ## methyHiC.snpwithmr.p.sort.vcf
out="k562"  ## methyHiC.methsnv.p.fragment_file
~/software/HapCUT2/build/extractHAIRS --hic 1 --bam ${input} --VCF ${vcf} --out ${out}.fragment_file
~/software/HapCUT2/build/HAPCUT2 --hic 1 --fragments ${out}.fragment_file --vcf ${vcf} --output ${out}.haplo.txt
#### other samples is same as above discribed.
### haplotype length
sh haplotype.len.sh
### asm analysis
~/software_devp/methyhaplo/bin/methyhaplo -M asm -a Y -m methratio.txt -s bs.align.sam -o mr.asmhap

## The distribution of ASM on TSS/TES/Exon start etc.
### The distribution of ASM on TSS
awk -v OFS="\t" '{if($7=="+"){print $1,$4,$7}else{print $1,$5,$7}}' ~/practice/Genome/hg38/gene.gtf > gene.ped3
~/software_devp/methyhaplo/ASMannoSites -o asmonTSS -G ~/practice/Genome/hg38/batmeth2-chr/hg38.chr.fa --ped ~/practice/Genome/hg38/gene.ped3 -d 5000 -ap newasm.plus.txt -an newasm.neg.txt
etc....
### The distribution of ASM on Exon start
awk -v OFS="\t" '$3=="exon"{if($7=="+"){print $1,$4,$7}else{print $1,$5,$7}}' ~/practice/Genome/hg38/Homo_sapiens.GRCh38.90.addchr.gtf > ~/practice/Genome/hg38/exon.ped3
~/software_devp/methyhaplo/ASMannoSites -o asmonExon -G ~/practice/Genome/hg38/batmeth2-chr/hg38.chr.fa -b ~/practice/Genome/hg38/exon.ped3 -d 5000 -ap newasm.plus.txt -an newasm.neg.txt
etc....
### The distribution of ASM on Exon end
awk -v OFS="\t" '$3=="exon"{if($7=="+"){print $1,$5,$7}else{print $1,$4,$7}}' ~/practice/Genome/hg38/Homo_sapiens.GRCh38.90.addchr.gtf > ~/practice/Genome/hg38/exonEnd.ped3
~/software_devp/methyhaplo/ASMannoSites -o asmonExonend -G ~/practice/Genome/hg38/batmeth2-chr/hg38.chr.fa -b ~/practice/Genome/hg38/exonEnd.ped3 -d 5000 -ap newasm.plus.txt -an newasm.neg.txt
etc...

## the distribution of ASM site in different functional components
### get hg38_intergenic without promoter (2Kb)
awk 'BEGIN{OFS="\t";}ARGIND==1 {a[$1]=$2} ARGIND==2 && $3=="gene" {if($7=="+"){if($4>2000){s=$4-2000}else{s=$4};print $1,s,$5}else{if($5+2000 <a[$1]){s=$5+2000}else{s=$5};print $1,$4,s}}' hg38.chr.fa.len.sort Homo_sapiens.GRCh38.90.addchr.gtf | sortBed | complementBed -i stdin -g hg38.chr.fa.len.sort > hg38_intergenic.noPromoter.bed
awk -v OFS="\t" '$1!~/^*/{print $1,$2"\n"$1,$3}' newasm.plus.txt  newasm.neg.txt | sort -k1,1 -k2,2n |uniq > asm.loci.txt
### intron
awk -v OFS="\t" '{print $1,$2,$2}' asm.loci.txt |bedtools intersect -a -  -b ~/practice/Genome/hg38/hg38_intron.bed -wo | awk '{a[$1" "$2]}END{print length(a)}'
### exon
awk -v OFS="\t" '{print $1,$2,$2}' asm.loci.txt |bedtools intersect -a -  -b ~/practice/Genome/hg38/hg38_exon_merged.bed -wo | awk '{a[$1" "$2]}END{print length(a)}'
### intergenic
awk -v OFS="\t" '{print $1,$2,$2}' asm.loci.txt |bedtools intersect -a -  -b ~/practice/Genome/hg38/hg38_intergenic.noPromoter.bed -wo | awk '{a[$1" "$2]}END{print length(a)}'
### Promoter
awk -v OFS="\t" '{print $1,$2,$2}' asm.loci.txt |bedtools intersect -a -  -b ~/practice/Genome/hg38/hg38_Promoter.bed -wo | awk '{a[$1" "$2]}END{print length(a)}'

## asm and gene expression
### alignment and get gene expression count
samtools sort -n RNAseq.sam -o RNAseq.sortname.bam 
htseq-count -f bam RNAseq.sortname.bam Homo_sapiens.GRCh38.90.addchr.gtf > RNAseq.Allgenecount
### divide to two cat
 awk 'ARGIND==1 && $2>0{a[$1]}ARGIND==2{id=$10;gsub(/"|;/, "", id);if(id in a){print $0}}' RNAseq.Allgenecount gene.gtf > gex.big0.gtf
$ awk 'ARGIND==1 && $2==0{a[$1]}ARGIND==2{id=$10;gsub(/"|;/, "", id);if(id in a){print $0}}' RNAseq.Allgenecount gene.gtf > gex.eq0.gtf
awk -v OFS="\t" '{if($7=="-"){print $1,$5,"-"}else{print $1,$4,"+"}}' gex.big0.gtf > gex.expression.TSS.ped3
awk -v OFS="\t" '{if($7=="-"){print $1,$5,"-"}else{print $1,$4,"+"}}' gex.eq0.gtf > gex.unexpression.TSS.ped3
### enrichment analysis
awk -v num=21238 '{printf $1;for(i=2;i<=NF;i++){printf "\t"100*$i/num}}' A549.asmonExpressionGene.Methy.1.txt > A549.asmonExpressionGene.Methy.1.txt.Aver
awk -v num=37005 '{printf $1;for(i=2;i<=NF;i++){printf "\t"100*$i/num}}' A549.asmonUnGene.Methy.1.txt > A549.asmonUnGene.Methy.1.txt.Aver


## ASM Accuracy Verification amrfinder/methpipe
### All chromosome
allelicmeth -c /public/home/qwzhou/project/Genome/arabidopsis/chrformethpipe/2.fa -o BS-cviler.allelicmeth.2 BS-cviler.epiread.2
allelicmeth -c /public/home/qwzhou/project/Genome/arabidopsis/chrformethpipe/3.fa -o BS-cviler.allelicmeth.3 BS-cviler.epiread.3
allelicmeth -c /public/home/qwzhou/project/Genome/arabidopsis/chrformethpipe/4.fa -o BS-cviler.allelicmeth.4 BS-cviler.epiread.4
allelicmeth -c /public/home/qwzhou/project/Genome/arabidopsis/chrformethpipe/5.fa -o BS-cviler.allelicmeth.5 BS-cviler.epiread.5
allelicmeth -c /public/home/qwzhou/project/Genome/arabidopsis/chrformethpipe/chloroplast.fa -o BS-cviler.allelicmeth.chloroplast BS-cviler.epiread.chloroplast
allelicmeth -c /public/home/qwzhou/project/Genome/arabidopsis/chrformethpipe/mitochondria.fa -o BS-cviler.allelicmeth.mitochondria BS-cviler.epiread.mitochondria
cat BS-cviler.allelicmeth.* > BS-cviler.allelicmeth
### methpipe
awk 'ARGIND==1 && $5<0.05 && sqrt($NF*$NF)>0.6' BS.cvi.ler.dmc | awk -v lp=0,lc="N",a="",an=0 '
{
    if(lc!=$1){
        if(an>1){
            print a;
        }
        lp=$2;lc=$1;a=$0;an=1
    }else{
        if(lp>=$2-300){
            a=a"\n"$0;an++
        }else if(an>1){
            print a;a="";
            a=$0;an=1;
        }else{a=$0;an=1;}
        lp=$2;lc=$1;
    }
}
' | awk 'ARGIND==1{a[$1" "$2]}ARGIND==2 && $5<0.0001{if($1" "($2+1) in a){print $0}}' - BS-cviler.allelicmeth |wc -l
7207
[17:10:00] qwzhou@c06n09:~/project/methhap/tair :
$ awk '$5<0.0001' BS-cviler.allelicmeth |wc -l
22741

### methyhaplo
~/software/methyhaplo/ASM -i BS-cviler.hap.cg.neg.txt -o BS-cviler.hap.asm.cg.neg.txt -p 0.05 -c 3
awk -v OFS="\t" '$1!~/^*/{print $1,$2"\n"$1,$3}' BS-cviler.asm.cg.plus.txt  | uniq > BS-cviler.asm.cg.plus.txt2
awk -v OFS="\t" '$1!~/^*/{print $1,$2"\n"$1,$3}' BS-cviler.asm.cg.neg.txt  | uniq > BS-cviler.asm.cg.neg.txt2
wc -l BS-cviler.asm.cg.plus.txt2
7306 BS-cviler.asm.cg.plus.txt2
wc -l BS-cviler.asm.cg.neg.txt2
6916 BS-cviler.asm.cg.neg.txt2
### plus
awk 'ARGIND==1 && $5<0.05 && sqrt($NF*$NF)>0.6' BS.cvi.ler.dmc | awk -v lp=0,lc="N",a="",an=0 '
{
    if(lc!=$1){
        if(an>1){
            print a;
        }
        lp=$2;lc=$1;a=$0;an=1
    }else{
        if(lp>=$2-300){
            a=a"\n"$0;an++
        }else if(an>1){
            print a;a="";
            a=$0;an=1;
        }else{a=$0;an=1;}
        lp=$2;lc=$1;
    }
}
' | awk 'ARGIND==1{a[$1" "$2]}ARGIND==2 && $1!~/^*/{if($1" "$2 in a){print $0}}' - BS-cviler.asm.cg.neg.txt2 |wc -l
4787
#### neg
awk 'ARGIND==1 && $5<0.05 && sqrt($NF*$NF)>0.6' BS.cvi.ler.dmc | awk -v lp=0,lc="N",a="",an=0 '
{
    if(lc!=$1){
        if(an>1){
            print a;
        }
        lp=$2;lc=$1;a=$0;an=1
    }else{
        if(lp>=$2-300){
            a=a"\n"$0;an++
        }else if(an>1){
            print a;a="";
            a=$0;an=1;
        }else{a=$0;an=1;}
        lp=$2;lc=$1;
    }
}
' | awk 'ARGIND==1{a[$1" "$2]}ARGIND==2 && $1!~/^*/{if($1" "$2 in a){print $0}}' - BS-cviler.asm.cg.plus.txt2 |wc -l
5137

--------------------------------------
## Statistical analysis of the overlap sites bettween wgbs hap and wgs hap in k562, hepG2
###snp hap
awk 'BEGIN{split("",s1);split("",d1);split("",s2);split("",d2);}
    ARGIND==1{a[$1" "$2]=$4;b[$1" "$2]=$5}ARGIND!=1{
    if($1~/^B/ && FNR>1){
        if(length(s1)>length(s2)){
            for(i in s1){same[i]=s1[i]};for(i in d1){diff[i]=d1[i]}
        }else{
            for(i in s2){same[i]=s2[i]};for(i in d2){diff[i]=d2[i]};
        }
        split("",s1);split("",d1);split("",s2);split("",d2);
    }else if($2==0 || $2==1){
        if($2==0){
            if($4" "$5 in a){
                if(a[$4" "$5]==$6){s1[$4" "$5]}else{if($8~/:/){d1[$4" "$5]}}
                if(b[$4" "$5]==$6){s2[$4" "$5]}else{if($8~/:/){d2[$4" "$5]}}
            }
        }else{
            if($4" "$5 in a){
                if(a[$4" "$5]==$7){s1[$4" "$5]}else{if($8~/:/){d1[$4" "$5]}}
                if(b[$4" "$5]==$7){s2[$4" "$5]}else{if($8~/:/){d2[$4" "$5]}}
            }
        }
    }
}END{
    if(length(s1)>length(s2)){
        for(i in s1){same[i]=s1[i]};for(i in d1){diff[i]=d1[i]}
    }else{
        for(i in s2){same[i]=s2[i]};for(i in d2){diff[i]=d2[i]};
    }
    for(i in same){print i,"same"};for(i in diff){print i,"diff"}
}' k562.wgshap.bed k562.snp.hapcut2.haplo.txt | sort -k1,1 -k2,2n > k562.snphap.wgshap.txt
### snpmr:
k562.wgshap.bed k562.snpwithmr.p.hapcut2.haplo.txt k562.snpwithmr.n.hapcut2.haplo.txt | sort -k1,1 -k2,2n > k562.mrsnphap.wgshap.txt
### merge file
awk -v OFS="\t" 'ARGIND==1{a[$1"\t"$2]="N\tN"}ARGIND==2{a[$1"\t"$2]=$3"\tN"}ARGIND==3{split(a[$1"\t"$2],b,"\t"); a[$1"\t"$2]=b[1]"\t"$3}END{for(i in a){print i,a[i]}}' k562.wgshap.bed k562.mrsnphap.wgshap.txt k562.snphap.wgshap.txt |sort -k1,1 -k2,2n > k562.vs.wgshap.merge.txt

### process file format
awk -v OFS="\t" '{print "chr"$1,$2,$2+1,$3,$4}' HepG2.haplotyped_across_tumor.passing.clean.txt > HepG2.haplotyped.clean.bed
### run liftover：
../liftOver HepG2.haplotyped.clean.bed hg19ToHg38.over.chain.gz HepG2.hapwgs.bed hepG2.unlifted.bed

### snp+meth
HepG2.hapwgs.bed hepG2.snp.hapcut2.haplo.txt > hepG2.snphap.wgshap.txt
HepG2.hapwgs.bed hepG2.snpwithmr.p.hapcut2.haplo.txt hepG2.snpwithmr.n.hapcut2.haplo.txt > hepG2.mrsnphap.wgshap.txt
### merge file
awk -v OFS="\t" 'ARGIND==1{a[$1"\t"$2]="N\tN"}ARGIND==2{a[$1"\t"$2]=$3"\tN"}ARGIND==3{split(a[$1"\t"$2],b,"\t"); a[$1"\t"$2]=b[1]"\t"$3}END{for(i in a){print i,a[i]}}' HepG2.hapwgs.bed hepG2.mrsnphap.wgshap.txt hepG2.snphap.wgshap.txt |sort -k1,1 -k2,2n > hepG2.vs.wgshap.merge.txt
awk '{print $1,$2,$3,$4}' hepG2.mrsnphap.wgshap.txt | awk '{if($3=="same"){a[$1" "$2]}else if(!($1" "$2 in a)){b[$1" "$2]}}END{print length(a),length(b)}'
sam/diff： 125816 3907
wgs hap snp: 370768，hetero-SNP：182255
snp in bs-seq: 1942278
snp in bs with wgs hap: 151546
snp hap with wgs hap: 80262, sam/diff:72551 7711
mrsnp hap with wgs hap: 129723，sam/diff：125816 3907

----------------------------------------------------------
## haplotype length in Arabidopsis thaliana
### merge plus and neg strand results
awk -v OFS="\t" '$1!~/\*/{if($1~/BLOCK/){if(chr!="" && end>start && start!=0){print "=== "chr,start,end"\n"inf};chr="";start=0;end=0;inf="";}else{if($2=="-"){if(chr!="" && end>start && start!=0){print "=== "chr,start,end"\n"inf};chr="";start=0;end=0;inf=""}else if($2==0 || $2==1){if(start==0){chr=$4;start=$5;inf=$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}else{inf=inf"\n"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7;if($5>end){end=$5}}}}}' BS-cviler.p.hapcut2.haplo.txt BS-cviler.n.hapcut2.haplo.txt  | awk -v OFS="\t" '$1~/^=/{print $2,$3,$4}' | bedtools merge | less
awk -v OFS="\t" 'sub(/Chr/,"",$1) && $4!~/transposon_fragment/' TAIR10_Transposable_Elements.bed > TAIR10_Transposable_Elements.replacechr.bed
### haplotype on repeat
#### 1. DNA meth + SNP
awk -v OFS="\t" '$1!~/\*/{if($1~/BLOCK/){if(chr!="" && end>start && start!=0){print "=== "chr,start,end"\n"inf};chr="";start=0;end=0;inf="";}else{if($2=="-"){if(chr!="" && end>start && start!=0){print "=== "chr,start,end"\n"inf};chr="";start=0;end=0;inf=""}else if($2==0 || $2==1){if(start==0){chr=$4;start=$5;inf=$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}else{inf=inf"\n"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7;if($5>end){end=$5}}}}}' BS-cviler.p.hapcut2.haplo.txt BS-cviler.n.hapcut2.haplo.txt  | awk -v OFS="\t" '$1~/^=/{print $2,$3,$4}' | sort -k1,1 -k2,2n -k3,3n | bedtools merge | bedtools intersect -a - -b TAIR10_Transposable_Elements.replacechr.bed -wo | awk '{a[$7]}END{print length(a)}'
#### 2. SNP
awk -v OFS="\t" '$1!~/\*/{if($1~/BLOCK/){if(chr!="" && end>start && start!=0){print "=== "chr,start,end"\n"inf};chr="";start=0;end=0;inf="";}else{if($2=="-"){if(chr!="" && end>start && start!=0){print "=== "chr,start,end"\n"inf};chr="";start=0;end=0;inf=""}else if($2==0 || $2==1){if(start==0){chr=$4;start=$5;inf=$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}else{inf=inf"\n"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7;if($5>end){end=$5}}}}}' BS-cviler.hapcut2.snp.haplo.txt  | awk -v OFS="\t" '$1~/^=/{print $2,$3,$4}' | sort -k1,1 -k2,2n -k3,3n | bedtools intersect -a - -b TAIR10_Transposable_Elements.replacechr.bed -wo | awk '{a[$7]}END{print length(a)}'
#### 3. DNA + SNP  (contain SNP)
awk -v OFS="\t" '$1!~/\*/{if($1~/BLOCK/){if(chr!="" && end>start && start!=0){print "=== "chr,start,end"\n"inf};chr="";start=0;end=0;inf="";}else{if($2=="-"){if(chr!="" && end>start && start!=0){print "=== "chr,start,end"\n"inf};chr="";start=0;end=0;inf=""}else if($2==0 || $2==1){if(start==0){chr=$4;start=$5;inf=$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}else{inf=inf"\n"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7;if($5>end){end=$5}}}}}' BS-cviler.p.hapcut2.haplo.txt BS-cviler.n.hapcut2.haplo.txt  | awk -v OFS="\t" '$1~/^=/{print $2,$3,$4}' | sort -k1,1 -k2,2n -k3,3n | bedtools merge | bedtools intersect -a - -b BS-cviler.snp.hetero.vcf -wo | awk -v OFS="\t" '{print $1,$2,$3}' |bedtools intersect -a - -b TAIR10_Transposable_Elements.replacechr.bed -wo | awk '{a[$7]}END{print length(a)}'

### ASM distribution on the genome
awk -v OFS="\t" 'gsub(/Chr/,"",$1) && $3!="chromosome"' ~/practice/Genome/arabidopsis/TAIR10_GFF3_genes.gff > ~/practice/Genome/arabidopsis/TAIR10_GFF3_genes.replacechr.gff
awk -v OFS="\t" '{if($7=="+"){if($4>2000){print $1,$4-2000,$4}else{print $1,1,$4}}else if($7=="-"){print $1,$5,$5+2000} }' ~/practice/Genome/arabidopsis/TAIR10.gene.modify.gff > ~/practice/Genome/arabidopsis/TAIR10.gene.modify.2Kpromoter.bed
cat BS-cviler.asm.plus.bed BS-cviler.asm.neg.bed | sort -k1,1 -k2,2n |bedtools merge | bedtools intersect -a - -b ~/practice/Genome/arabidopsis/TAIR10_GFF3_genes.replacechr.gff -f 0.50 -wo | awk '{a[$6" "$7" "$8]}END{for(i in a){print i}}' | awk '{a[$1]++}END{for(i in a){print i,a[i]}}'
ncRNA 54
pseudogene 69
exon 4174
rRNA 1
protein 3195
tRNA 3
transposable_element_gene 898
mRNA_TE_gene 900
five_prime_UTR 75
three_prime_UTR 133
pseudogenic_transcript 69
mRNA 3776
pseudogenic_exon 61
CDS 2932
gene 3041
Promoter 1890
TE 1588
cat BS-cviler.asm.plus.bed BS-cviler.asm.neg.bed | sort -k1,1 -k2,2n |bedtools merge | bedtools intersect -a - -b ~/practice/Genome/arabidopsis/TAIR10.gene.modify.2Kpromoter.bed -f 0.50 -wo | awk '{a[$4" "$5" "$6]}END{print length(a)}'
1890
cat BS-cviler.asm.plus.bed BS-cviler.asm.neg.bed | sort -k1,1 -k2,2n |bedtools merge | bedtools intersect -a - -b ~/practice/Genome/arabidopsis/TAIR10_Transposable_Elements.replacechr.bed -f 0.50 -wo | awk '{a[$4" "$5" "$6]}END{print length(a)}'
1588
### ASM distribution across the gene
~/software_devp/methyhaplo/ASManno -o ASMonGene -G ~/practice/Genome/arabidopsis/arabidopsis_batmeth2_index/TAIR10_chr_all.fa -gff ~/practice/Genome/arabidopsis/TAIR10.protein_coding_gene.gff -ap BS-cviler.asm.plus.bed -an BS-cviler.asm.neg.bed
~/practice/A549_BS/WGBS-HAIR/ASManno -o ASMonGene -G ~/practice/Genome/arabidopsis/arabidopsis_batmeth2_index/TAIR10_chr_all.fa -gff ~/practice/Genome/arabidopsis/TAIR10.protein_coding_gene.gff -ap BS-cviler.asm.plus.bed -an BS-cviler.asm.neg.bed -s 0.01