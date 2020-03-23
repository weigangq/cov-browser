#!/bin/bash

## GENERATE HAPLOTYPES ##
## output: A vcf file and a fasta sequence file of haplotypes ##
human_file=$1
ref_seq=$2
ploidy=$3
environ_files=$4

if [ $# -ne 4  ]
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

## Get sample FASTA
bcftools query -l snps2.vcf > samples
cat samples | while read line; do echo ">$line"; bcftools query -s "$line" -f '[%TGT]' snps2.vcf; echo; done > samples.fas

## check the number of samples in fas file and vcf file
bioseq -n samples.fas >> human.log
bcftools stats snps2.vcf > stats 
grep "number of samples" stats >> human.log

#vcftools --vcf --counts snps.vcf

## to get ref FASTA; if necessary
#echo ">ref" >> samples.fas
#grep "^CP" gbs50-snps2.vcf | cut -f4 | paste -s -d '' >> sample.fas)

#rename sequence
#cat samples.fas | sed "s/EPI_IS//; s/.sorted.bam//" > samples2.fas
#cat samples.fas | sed "s/.sorted.bam//" > samples2.fas
cat samples.fas | sed "s|^[^/]*/|>|g; s/.sorted.bam//" > samples2.fas

#get max haplotype 
max=$(bioseq -g samples.fas | bioseq -l | sort -k2 -nr | less | head -1 | cut -f2)
cutoff=$(echo $(( max*20/100 )))

## remove seqs with 20% of missing/non-ATCG bases "."
#bioseq -d'ambig:20' samples2.fas > samples3.fas
if [ $cutoff -gt 1 ]
then
	bioseq -d "ambig:$cutoff" samples2.fas > samples3.fas
else 
	bioseq -d "ambig:1" samples2.fas > samples3.fas
fi

#cat samples2.fas | grep "N" | grep -v

bioaln -i'fasta' --gap-char '-' samples3.fas  > samples4.aln

# change gap to "-" (internal gaps counted as different)
#bioaln -i'fasta' --gap-char '-' samples3.fas  > samples4.aln

# change gaps to "n" to help reducing unique haps
bioaln -i'fasta' --gap-char 'n' samples3.fas  > samples4.aln

