## Create alignment pileup and call variants using plodity file(plodity 1), multiallelic, first output is bam then piped to bcf ##
## change the directory of the sorted.bam files according for your run ##
# or -P 0.1; large P value for less strict call, default 1.1e-3 ##
bcftools mpileup -Ou -f ref.fas  /home/sakther/Borrelia_Thesis_Projects/corona_virus_study/cov-may-22-2020/human/*sorted.bam  /home/weigang/cov-may-22-2020/human/*sorted.bam /home/lli/covid/human/*sorted.bam | bcftools call -mv --ploidy-file ploidy.txt  -Ob -o calls.bcf -P 0.1

## Get Stats and check TS/TV ratio ( expeced more transitions than transversions) #for cov2 should be around more or less 2.5 
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
bcftools query -l snps2-ref.vcf > samples
sed -E "s/.+(EPI_ISL_[0-9]+).+.sorted.bam/\1/" samples > samples-short.txt
#change the name of the sample (verify this part since I did bit differently in my run)
bcftools reheader -s samples-short.txt snps2.vcf > snps3.vcf
cat samples-short.txt | while read line; do echo ">$line"; bcftools query -s "$line" -f '[%TGT]' snps3.vcf; echo; done > samples.fas

## check the number of samples in fas file and vcf file
bioseq -n samples.fas >> human.log
bcftools stats snps3.vcf > stats
grep "number of samples" stats >> human.log

#run the impute-hap.pl to count the # of missing bases in each SNP position 
#change the location of the impute-hap.pl file accordingly
perl /home/sakther/cov-browser/scripts/impute-hap.pl --dump-missing --format 'fasta' samples.fas > missing
sort -n missing | uniq -c | head #note the SNP position at the begining of the alignment from where we started to see very few missing bases  
sort -n missing | uniq -c | tail #note the SNP position at the end of the alignment from where we see very few missing bases 

#covert the gap to "n"
bioaln -i'fasta' --gap-char 'n' samples.fas > samples-for-no-impute.aln
#check the # of sequences 
bioaln -n samples-for-no-impute.aln

#change the location of the impute-hap.pl file accordingly
#for my case I removed the positions before 6 and after 125 (those positions have way too many missing bases)
perl /home/sakther/cov-browser/scripts/impute-hap.pl --start 6 --end 125 --format 'clustalw' --no-imputation samples-for-no-impute.aln 2> no-impute.log

# the above command output 3 files 
# a log file which we will use to map the sample ids with STs
# a alignment file of unique STs
# a alignment file of all the samples 

#check the number of sequences and SNPs in both alignment files 
bioaln -n no-imputed-sample.aln #you can use bioseq -l as well 
bioaln -n no-imputed.aln

#remove the "filtered position of the alignment" from the vcf file
# find the absolute position of filtered start and end in the vcf file 
grep "^NC" snps3.vcf | cut -f2 | cat -n | head #in my case I looked for the position of the SNP #6 in vcf file which is pos 187
grep "^NC" snps3.vcf | cut -f2 | cat -n | tail #in my case I looked for the position of the SNP #125 in vcf file which is pos 29734

vcftools --chr NC_045512 --from-bp 187 --to-bp 29734  --recode --recode-INFO-all --vcf snps3.vcf --out snps4
mv snps4.recode.vcf snps4.vcf

#validate the length of the alignment file and new vcf file. Both should be same 
#check the # of SNPs in stats
bcftools stats snps4.vcf > stats
bioaln -l no-imputed.aln













