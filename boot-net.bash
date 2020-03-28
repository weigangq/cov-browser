#!/usr/bin/env bash
# Usage: $0 <boot number>
dir='boot-files';
end=$1;
# bootstrap alignment, build MST, and get edges
perl hapnet.pl --genome ref.gb --vcf snps2-ref-03-19.vcf --hap imputed.aln --impute-log impute.log --output edge > $dir/net-0.edge
for i in $(seq 1 $end); do 
    echo -ne "boot $i\t";
    bioaln --boot imputed.aln > $dir/boot-$i.aln; 
    perl hapnet.pl --genome ref.gb --vcf snps2-ref-03-19.vcf --hap $dir/boot-$i.aln --impute-log impute.log --output edge > $dir/net-$i.edge
    cut -f1,2 $dir/net-$i.edge | sort -u | tr '\t' '-' > $dir/edge-$i; 
    echo -ne "done boot $i\n";
done

# obtain majority-rule parents
cat $dir/edge-* | sort | uniq -c | sort -rn > unique-edges
perl majority-parent.pl unique-edges > pa-cts
exit;
