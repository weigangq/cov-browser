#!/bin/bash

## GENERATE HAPLOTYPES ##
## output: A vcf file and a fasta sequence file of haplotypes ##

human_file=$1 #a fasta file contain all the human cov genome sequences
ref_seq=$2
ploidy=$3 #vcf ploidy file
environ_files=$4 #folder name that contain non-human cov genome bam files
out_file=$5 #outgroup bam file (compressed)

if [ $# -ne 5  ]
  then
    echo "Not enough arguments provided"
	exit 1
fi

#remove files and directory for previous run (if any) 
if [ -e human.log ];then rm human.log ; fi
if [ -e human ];then rm -rf human ; fi 

##funcion sam align generates alignment with nucmer, sam, and indexed sorted bam##
sam_align() {
fas=$1;
ref=$2;
bin_path=/home/weigang/mummer-4.0.0beta2;
#ref_path=/home/weigang/cov-03-09-2020;
echo "processing $fas";
id=$(basename $fas .fas);
$bin_path/nucmer --sam-long=$id.sam $ref $fas;
samtools view -b $id.sam -T $ref > $id.bam;
samtools sort $id.bam -o $id.sorted.bam;
samtools index $id.sorted.bam;
if [ ! -s "$id.sorted.bam" ] 
then
	echo $id.sorted.bam >> human.log
fi
echo "done $fas"
}

#make sure you don't have a directory of same name. Otherwise it will get deleted  
mkdir human;
cd human;

## Break sequence file into individual files##
bioseq -B ../$human_file

## call funcion sam align to generate alignment with nucmer, sam, and indexed sorted bam ##

for f in *.fas
	do sam_align $f ../$ref_seq
done;

#back to original directory
cd ..

## Create alignment pileup and call variants using plodity file(plodity 1), multiallelic, first output is bam then piped to bcf 
bcftools mpileup -Ou -f $ref_seq   human/*.sorted.bam  $environ_files/*sorted.bam | bcftools call -mv --ploidy-file ploidy.txt  -Ob -o calls.bcf -P 0.05 ## or -P 0.1; large P value for less strict call, default 1.1e-3

## Get Stats and check TS/TV ratio ( expeced more transitions than transversions)
bcftools stats calls.bcf > stats
grep "^TSTV" stats >> human.log

##  get only biallelic SNPs ##
bcftools view -m2 -M2 --types snps calls.bcf > snps.bcf

## to verify run the stats again 
bcftools stats snps.bcf > stats
grep "^TSTV" stats >> human.log

## bcf to vcf
bcftools view snps.bcf > snps.vcf

## filter sites by allele counts: only keeps informative sites
vcftools --vcf snps.vcf --mac 2 --recode --recode-INFO-all --out snps2.vcf

## rename the snp2 recode file 
mv snps2.vcf.recode.vcf snps2.vcf

#step 1: index with gzip & tabix
if [ -e snps2.vcf.gz ];then rm snps2.vcf.gz ; fi
bgzip snps2.vcf
tabix snps2.vcf.gz

##fitler outgroup SNPs
bcftools isec -n=2 -w1 $out_file snps2.vcf.gz > otg-isec-sites.vcf

## index with gzip & tabix
if [ -e otg-isec-sites.vcf.gz ];then rm otg-isec-sites.vcf.gz ; fi
bgzip otg-isec-sites.vcf
tabix otg-isec-sites.vcf.gz

## merge (no duplicated line; add ref state "0" instead of missing)
vcf-merge -d -R 0 otg-isec-sites.vcf.gz snps2.vcf.gz > snps-ref.vcf

## check stats 
bcftools stats snps-ref.vcf > stats
grep "number of samples" stats >> human.log

##rename sequecnes 
sed -e 's/\(human\/\|environ\/\|.sorted.bam\)//g' snps-ref.vcf > snps2-ref.vcf

## Get sample FASTA
bcftools query -l snps2-ref.vcf > samples
cat samples | while read line; do echo ">$line"; bcftools query -s "$line" -f '[%TGT]' snps2-ref.vcf; echo; done > samples.fas

## check the number of samples in fas file and vcf file
bioseq -n samples.fas >> human.log
bcftools stats snps2-ref.vcf > stats 
grep "number of samples" stats >> human.log

#get max haplotype 
max=$(bioseq -g samples.fas | bioseq -l | sort -k2 -nr | less | head -1 | cut -f2)
cutoff=$(echo $(( max*10/100 )))
echo "ambig cutoff is $cutoff" >> human.log

## remove seqs with 20% of missing/non-ATCG bases "."
if [ $cutoff -gt 1 ]
then
	bioseq -d "ambig:$cutoff" samples.fas > samples2.fas
else 
	bioseq -d "ambig:1" samples.fas > samples2.fas
fi

# change gap to "-" (internal gaps counted as different)
#bioaln -i'fasta' --gap-char '-' samples3.fas  > samples4.aln

# change gaps to "n" to help reducing unique haps
bioaln -i'fasta' --gap-char 'n' samples2.fas  > samples3.aln
