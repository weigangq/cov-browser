#!/usr/bin/env bash
# Usage: $0 <boot number>
dir_file=$HOME/boot-files;
dir_exe=$HOME/cov-browser;
end=$1;

rm $dir_file/*;
# perl hapnet.pl --genome ref.gb --vcf snps-ref.vcf --impute-log impute.log --hap imputed.aln
#1 get edges for the original alignment ("imputed.aln")
perl $dir_exe/hapnet.pl --genome $dir_exe/ref.gb --vcf snps-ref.vcf --hap imputed2.aln --impute-log impute.log2 2> log
mv net.json $dir_file/net-0.json;
mv edges.tsv $dir_file/edges-0.tsv;
mv nodes.tsv $dir_file/nodes-0.tsv;
cut -f1,2 $dir_file/edges-0.tsv | sort -u | tr '\t' '-' > $dir_file/edge-0;

#2 get edges for bootstraped alignments
for i in $(seq 1 $end); do 
    echo -ne "boot $i\t";
    bioaln --boot imputed2.aln > $dir_file/boot-$i.aln; 
    perl $dir_exe/hapnet.pl --genome $dir_exe/ref.gb --vcf snps-ref.vcf --hap $dir_file/boot-$i.aln --impute-log impute.log2 2> log
    mv net.json $dir_file/net-$i.json;
    mv edges.tsv $dir_file/edges-$i.tsv;
    mv nodes.tsv $dir_file/nodes-$i.tsv;
    cut -f1,2 $dir_file/edges-$i.tsv | sort -u | tr '\t' '-' > $dir_file/edge-$i; 
    echo -ne "done boot $i\n";
done

#3 obtain majority-rule parents
cat $dir_file/edge-* | sort | uniq -c | sort -rn > unique-edges
perl $dir_exe/majority-parent.pl unique-edges > pa-cts.tsv
exit;
