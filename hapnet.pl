#!/usr/bin/env perl
# 3/15/2020
# calculate minimum spanning tree of a set of haplotypes
# Inputs: haplotype seq alignment in FASTA; VCF file; BED file (embeded)
# Output: JSON for D3-force rendering
# Done (3/15/2020):  add additional attributes including gene, codon, AA, syn/nonsyn (by parsing VCF & BED file)
# Done (3/15/2020):  uniq_seq: consolidate sequences with internal missing bases; solved bioaln --gap-char "n" (not ideal)
# TO DO:
#     outgroup rooting: re-orient MSP path (Single-Source Shortest Paths, SSSP, won't work because it doesn't cover all nodes)
#
# Perl modules for Graphs: https://metacpan.org/release/Graph
#

use strict;
use Bio::AlignIO;
use Data::Dumper;
use Graph;
use Graph::Writer::Dot;
use Graph::Traversal::DFS;
use Graph::Traverse;
use Bio::SeqFeature::Generic;
use JSON;
use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO;
use Bio::Tools::CodonTable;
my $myCodonTable   = Bio::Tools::CodonTable->new();
my $writer = Graph::Writer::Dot->new();
#die "Usage: $0 --hap <fasta> --vcf <vcf> --genome <fasta>\n" unless;
my %options;
my $outgroup = 'EPI_ISL_402131';
GetOptions (
    \%options,
    "hap=s",
    "vcf=s",
    "genome=s",
    "output=s",
    "help",
    ) or pod2usage(2);

pod2usage(1) if $options{'help'};
#print Dumper(\%options); exit;
#####################################
# Read Files
#######################################
# Read FASTA alignment
my $in = Bio::AlignIO->new(-file=>$options{'hap'}, -format=>'clustalw');
my $aln = $in->next_aln();
my $alnLen = $aln->length();
my ($stAln, $refSt) = $aln->uniq_seq(); # this doesn't count gaps on edges but still treats internal gaps as distinct seqs
my %STs;
foreach my $st ($stAln->each_seq) {
    $STs{$st->display_id()} = $st;
}

warn "Number of STs:", scalar keys %STs, "\n"; 

my %hapST;
foreach my $st (keys %$refSt) {
    $hapST{'ST'. $st} = $refSt->{$st}
}

#print Dumper(\%hapST); exit;
# Read VCF file
my @pos;
open VCF, "<", $options{'vcf'};
while(<VCF>){
    chomp;
    next unless /^NC_/;
    my @a = split;
    push @pos, $a[1];
}
close VCF;

my $gIN = Bio::SeqIO->new(-file => $options{'genome'}, -format => 'genbank');
my $ge = $gIN->next_seq();

# Read BED data & create genome features
my @gene_feats;
while (<DATA>) {
    chomp;
    my ($chr, $start, $end, $gid, $strand, $rest) = split;
    push @gene_feats, Bio::SeqFeature::Generic->new(-seq_id => $gid, -start => $start, -end => $end, -strand => $strand);
}

#print Dumper(\@gene_feats); 
my @igs_feats;
for (my $i=0; $i<$#gene_feats; $i++) {
    my $this = $gene_feats[$i];
    my $next = $gene_feats[$i+1];
    push @igs_feats, Bio::SeqFeature::Generic->new(-seq_id => "igs:" . $this->seq_id() . "|" . $next->seq_id(),
						   -start => $this->end() + 1, 
						   -end => $next->start() - 1, 
						   -strand => $this->strand() . "|" . $next->strand(),
	);    
}
#print Dumper(\@igs_feats); 
my @utr_feats;
push @utr_feats, Bio::SeqFeature::Generic->new(-seq_id => '5-UTR', 
					       -start => 1, 
					       -end => $gene_feats[0]->start -1, 
					       -strand => 1);

push @utr_feats, Bio::SeqFeature::Generic->new(-seq_id => '3-UTR', 
					       -start => $gene_feats[-1]->end + 1, 
					       -end => 29903, # end of NC_045512.2 
					       -strand => 1);

#print Dumper(\@utr_feats); exit;
#############################
# Make a graph & calculate MST & polorize with root
#############################
my $graph = Graph->new(undirected => 1);
my @sts = sort keys %STs;
for (my $i=0; $i<$#sts; $i++) { # build a complete graph from edges
    my $idI = $STs{$sts[$i]}->display_id();
    for (my $j=$i+1; $j<=$#sts; $j++) {
	my $idJ = $STs{$sts[$j]}->display_id();
	$graph->add_edge($idI, $idJ);
	my @diffs = &comp_seq($STs{$sts[$i]}, $STs{$sts[$j]});
#	print Dumper(\@diffs);
	$graph->set_edge_weight($idI, $idJ, scalar @diffs); # seq diff as weight (to be minimized in MST)
    }
}

my (@nodelist, @edgelist);
my $mstg = $graph->MST_Kruskal();  # Krustal's MST Algorithm
my @Vs = $mstg->unique_vertices();
my $root;

# Call back function to process each vertex using an anonymous sub
my $find_root_and_attach_haps = sub {
    my ($v, $self) = @_;
    my $ref = $hapST{$v};
    foreach (@$ref) { 
	if ($_ eq $outgroup) {
	    $root = $v;
	    last; 
	}
    }
    $mstg->set_vertex_attribute($v, 'haps', $hapST{$v});
    $mstg->set_vertex_attribute($v, 'id', $v);
    push @nodelist, $mstg->get_vertex_attributes($v);
};

# 1st traversal
my $traversal = Graph::Traversal::DFS->new(
    $mstg,
    pre => $find_root_and_attach_haps # process each vertex by pre-order traversal
    );
$traversal->dfs();

# Call back function to process each edge using an anonymous sub
my $polarize_an_edge = sub {
    my ($u, $v, $self) = @_;
    my @diffs = &comp_seq($STs{$u}, $STs{$v});
    $mstg->set_edge_attribute($u, $v, 'change', \@diffs);
    $mstg->set_edge_attribute($u, $v, 'source', $u);
    $mstg->set_edge_attribute($u, $v, 'target', $v);
    $mstg->set_edge_attribute($u, $v, 'id', join "-", ($u, $v));
    push @edgelist, $mstg->get_edge_attributes($u, $v);
};

# 2nd traversal
$traversal = Graph::Traversal::DFS->new(
    $mstg,
    start => $root,
    pre_edge => $polarize_an_edge
    );
$traversal->dfs();

#############################
# write out JSON
##############################

print to_json([\@nodelist, \@edgelist]) if $options{'output'} eq 'json';
$writer->write_graph($mstg, 'mst.dot') if $options{'output'} eq 'dot';
exit;

#####################################
# subroutines
###################################

sub comp_seq { # possibly no change? due to missing base? yes; fix uniq_seq
    my ($i, $j) = @_;
    my @str1 = split '', $i->seq();
    my @str2 = split '', $j->seq();
    my $id1 = $i->display_id();
    my $id2 = $j->display_id();
    my @dif;
    for (my $k=0; $k<$alnLen; $k++) {
	next unless $str1[$k] =~ /[atcg]/i; # i is missing, skip
	next unless $str2[$k] =~ /[atcg]/i; # j is missing, skip
	next if $str1[$k] eq $str2[$k]; # no diff, skip
	my $ref = &genome_info($pos[$k], $i->subseq($k+1,$k+1), $j->subseq($k+1,$k+1));
	push @dif, {'site' => $pos[$k], 
		    'src_base' => $i->subseq($k+1,$k+1), 
		    'tgt_base' => $j->subseq($k+1,$k+1),
		    'label' => $ref->[0],
		    'length' => $ref->[1],
		    'pos' => $ref->[2],
		    'src_codon' => $ref->[3] || undef,
		    'tgt_codon' => $ref->[4] || undef,
		    'src_aa' => $ref->[5] || undef,
		    'tgt_aa' => $ref->[6] || undef,
	}; # found a diff
    } 
    return @dif;
}

sub genome_info {
    my ($site, $base1, $base2) = @_;
    my $feat = &genome_region($site);
    if ($feat->seq_id =~ /UTR/ || $feat->seq_id =~ /igs/) { # IGS or UTR site
	return [ $feat->seq_id,  $feat->length, $site - $feat->start + 1 ];
    } else { # ORF site
	my $dna = Bio::Seq->new(-id=>'frag', -seq=>$ge->subseq($feat->start, $site));
	my $dnaLen = $dna->length();
	my $mod = $dnaLen % 3;
	my ($codon1, $codon2, $aa1, $aa2);
	if ($mod == 1) { # 1st codon position
	    $codon1 = $base1 . $ge->subseq($site+1, $site+2); 
	    $codon2 = $base2 . $ge->subseq($site+1, $site+2);
	} elsif ($mod == 2) {
	    $codon1 = $ge->subseq($site-1, $site-1) . $base1 . $ge->subseq($site+1, $site+1); 
	    $codon2 = $ge->subseq($site-1, $site-1) . $base2 . $ge->subseq($site+1, $site+1); 
	} else { # 3rd codon position
	    $codon1 = $ge->subseq($site-2, $site-1) . $base1; 
	    $codon2 = $ge->subseq($site-2, $site-1) . $base2; 
	}
	$aa1 = $myCodonTable->translate($codon1);
	$aa2 = $myCodonTable->translate($codon2);
	return [ $feat->seq_id,  $feat->length, $site - $feat->start + 1, $codon1, $codon2, $aa1, $aa2 ];
    }
}

sub genome_region {
    my $site = shift;
#    print $site, "=>";
    my $feat;
    foreach my $fe (@gene_feats, @igs_feats, @utr_feats) {
	next if $site < $fe->start;
	next if $site > $fe->end;
	$feat = $fe;
    }
#    print $feat->seq_id, ":", $feat->start, ":", $feat->end, "\n"; 
    return $feat;
}


__DATA__
NC_045512.2     265     21555   orf1ab  0       +
NC_045512.2     21562   25384   S       0       +
NC_045512.2     25392   26220   ORF3a   0       +
NC_045512.2     26244   26472   E       0       +
NC_045512.2     26522   27191   M       0       +
NC_045512.2     27201   27387   ORF6    0       +
NC_045512.2     27393   27759   ORF7a   0       +
NC_045512.2     27755   27887   ORF7b   0       +
NC_045512.2     27893   28259   ORF8    0       +
NC_045512.2     28273   29533   N       0       +
NC_045512.2     29557   29674   ORF10   0       +
