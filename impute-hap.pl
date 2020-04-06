#!/usr/bin/env perl
# 3/26/2020
# impute missing NTs based on closest nabe
#
# input: haplotypes from sam-align pipeline, each seq with < 10% gaps ("n", ".", or "-"; any non-ATCG)
#
# outputs: 
#   1. log file: removed seqs with 10% gaps & imputed haps
#   2. ST alignment
#   3. xref between STs and samples
#   4. retained sites
#
# Algorthm
#  1. remove seqs with 10% or more missing
#  2. remove terminal sites: 1-5 & 126-132 (need to determine by tabulate)
#  3. impuate middle missing bases ("n") by the nearest non-ATCG nabe
#  4. output unique STs

use strict;
use Bio::AlignIO;
use Data::Dumper;
use Getopt::Long;
use Bio::SimpleAlign;

my (%seqs, %hapSTs, %options);
GetOptions (
    \%options,
    "help",
    "debug|d",
    "no-imputation|n",
    "start=i",
    "end=i",
    "format=s",
    "dump-missing",
    ) or pod2usage(2);

pod2usage(1) if $options{'help'};
#die "Usage: $0 --start <n> --end <n> [--dump-missing] [--no-imputation] samples.aln\n" unless $options{'dump-missing'};
my $perc_missing = 0.1;
my $format = $options{format} || 'fasta';
my $in = Bio::AlignIO->new(-file => shift @ARGV, -format => $format);
my $aln = $in->next_aln();
my $len = $aln->length();
my @edges = ($options{start}||1, $options{end}||$len); # retain coords excludes between # to be determined by running "dump-missing | sort -n | uniq -c"
foreach my $seq ($aln->each_seq) {
    my @pos = @{&count_nonATCG($seq->seq)}; 

    if (scalar @pos >=  ($edges[1] - $edges[0] + 1) * $perc_missing) { # dump seqs with 10% or more non-ATCGs;
	warn "removed ", $seq->display_id, " for having >= 10% non-ATCGs: ", scalar @pos, "\n";
	next;
    } 

    $seqs{$seq->id()} = { 'seq' => $seq->trunc($edges[0], $edges[1]), 
			  'seq_ori' => $seq, 
			  'missing_pos' => \@pos,
			  'imputed_seq' => undef, # don't overwrite original seq
			  'imputed' => 0,
    }
}

&dump_missing if $options{'dump-missing'};

my @ids = sort keys %seqs;
#print Dumper(\@ids); exit;
#my %nabes;
# pairwise diff is the most costly with O(n2); don't scale well. Simplify by stop at the closest non-ATCG nabes?
my $ct = 0;
foreach my $id1 (@ids) {
    $ct++;
    next unless @{$seqs{$id1}->{'missing_pos'}}; # no ambig; no need to impute
    next if $options{'no-imputation'}; # do not impute
    my %nabes;
    foreach my $id2 (@ids) {
	next if @{$seqs{$id2}->{'missing_pos'}}; # has ambig; can't be used to impute
	my $pair = Bio::SimpleAlign->new();
        $pair->add_seq($seqs{$id1}->{seq}); # truncated
        $pair->add_seq($seqs{$id2}->{seq});
	my $diff = 100 - $pair->percentage_identity();
	$nabes{$id2} = $diff;
    }
    my @nabes = sort { $nabes{$a} <=> $nabes{$b} } keys %nabes;
    print STDERR "processing seq $ct ... \t";
    if (@nabes) {
	$seqs{$id1}->{'imputed'} = 1; 
	my $id2 = shift @nabes;
	$seqs{$id1}->{'imputed_seq'} = &impute_by_pos($seqs{$id1}, $seqs{$id2});
	warn "impute $id1 by ", $seqs{$id1}->{'imputed_seq'}->id(), "\n";
    } else {
	warn "not imputed: $id1; no nearest non-ATCG nabes\n" unless $seqs{$id1}->{'imputed'}
    }
}


########## Do not add seq obj; create a new obj!!!################################################
#	$seqs{$id}->{'imputed_seq'} = Bio::LocatableSeq->new(-seq => $seqs{$id2}->{seq}->seq,
#							     -id  => "$id|by-$id2"); # add seq
#################################################################################################

my $new = Bio::SimpleAlign->new();
foreach my $id (@ids) {
#    print $id, "\n";
#    my $seq = $seqs{$id}->{'seq'};
#    print $id, "\t", $seqs{$id}->{seq}->id(), "\t", $seqs{$id}->{'imputed'} ? $seqs{$id}->{'imputed_seq'}->id() : 'na', "\n";
    my $seq = $seqs{$id}->{'imputed'} ? $seqs{$id}->{'imputed_seq'} : $seqs{$id}->{'seq'};
    $new->add_seq($seq);
}
#$new->set_displayname_flat();
#my $out1 = Bio::AlignIO->new(-file=> ">new.aln", -format=>'clustalw');
#$out1->write_aln($new);

my ($st_aln, $ref_hap) = $new->uniq_seq();
$st_aln->set_displayname_flat();
#################
# output ST alignment & a log file (with sites and hap xrefs)
#################

my $out = Bio::AlignIO->new(-file=> ">imputed.aln", -format=>'clustalw');
$out->write_aln($st_aln);

open LOG, ">", "impute.log";
print LOG $edges[0], "\t", $edges[1], "\n";

my %haps = %$ref_hap;
foreach my $st (sort {$a <=> $b} keys %haps) {
    my @hap = @{$haps{$st}};
    foreach my $hap (@hap) {
	print LOG "ST", $st, "\t", $hap eq 'cov-outgroup' ? 'EPI_ISL_402131' : $hap;
	if ($hap =~ /(EPI_\S+)\|by-EPI_\S+/) {
	    my $id = $1;
	    print LOG "\t";
	    print LOG join("-", @{$seqs{$id}->{'missing_pos'}});
	    print LOG "\n";
	} else {
	    print LOG "\n";
	}
    }
}
close LOG;

exit;

sub impute_by_pos {
    my ($seq_1, $seq_2) = @_;
    my $seq_ob1 = $seq_1->{'seq_ori'};
    my $seq_ob2 = $seq_2->{'seq_ori'};
    my $seq_str1 = $seq_ob1->seq();
    my $seq_str2 = $seq_ob2->seq();
    my $id1 = $seq_ob1->id();
    my $id2 = $seq_ob2->id();
    my @pos = @{$seq_1->{'missing_pos'}};
    foreach my $p (@pos) {
	my $base1 = $seq_ob1->subseq($p, $p);
	my $base2 = $seq_ob2->subseq($p, $p);
	warn "$p: $base1 by $base2\n";
	substr($seq_str1, $p-1, 1, $base2);
    }
    my $seq_replace = Bio::LocatableSeq->new(-id => "$id1|by-$id2",
					     -seq => $seq_str1
	);

    return $seq_replace->trunc($edges[0], $edges[1]);
}


sub dump_missing {
    foreach my $id (keys %seqs) {
	next unless @{$seqs{$id}->{'missing_pos'}};
	print join "\n", @{$seqs{$id}->{'missing_pos'}};
	print "\n";
    }
    exit;
}


sub count_nonATCG {
    my $seq_str = shift;
    my @bases = split '',  $seq_str;
    my @miss;
    for(my $i=0; $i<=$#bases; $i++) {
	next if ($i <= $edges[0] - 2 || $i >= $edges[1]) && !$options{'dump-missing'};
	push @miss, $i+1 unless $bases[$i] =~ /[ATCG]/i;
    }
    return \@miss;
}
