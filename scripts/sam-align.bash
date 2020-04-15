#!/bin/bash
fas=$1;
bin_path=$HOME/mummer-4.0.0beta2;
ref_path=$HOME/cov-browser/ref-genome;
echo "processing $fas";
id=$(basename $fas .fas);
$bin_path/nucmer --sam-long=$id.sam $ref_path/ref.fas $fas;
samtools view -b $id.sam -T $ref_path/ref.fas > $id.bam;
samtools sort $id.bam -o $id.sorted.bam;
samtools index $id.sorted.bam;
echo "done $fas"
exit;
