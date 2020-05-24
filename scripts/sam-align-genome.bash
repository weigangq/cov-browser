#!/bin/bash

## GENERATE HAPLOTYPES ##
## output: A vcf file and a fasta sequence file of haplotypes ##

human_file=$1 #a fasta file contain all the human cov genome sequences
ref_seq=$2

if [ $# -ne 2  ]
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
#bin_path=/home/weigang/mummer-4.0.0beta2;
#ref_path=/home/weigang/cov-03-09-2020;
echo "processing $fas";
id=$(basename $fas .fasta);
nucmer --sam-long=$id.sam $ref $fas;
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

