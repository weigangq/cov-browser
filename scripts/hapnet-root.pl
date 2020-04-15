#!/usr/bin/env perl
# 3/15/2020
# calculate minimum spanning tree of a set of haplotypes
# Inputs: haplotype seq alignment in FASTA; VCF file; BED file (embeded)
# Output: JSON for D3-force rendering
# Done (3/15/2020):  add additional attributes including gene, codon, AA, syn/nonsyn (by parsing VCF & BED file)
# Done (3/15/2020):  uniq_seq: consolidate sequences with internal missing bases; solved bioaln --gap-char "n" (not ideal)
# Done (3/16/2020):  outgroup rooting: re-orient MSP path (Single-Source Shortest Paths, SSSP, won't work because it doesn't cover all nodes)
# To Do: Impute missing NTs? currently STs are represented by seqs encountered first (input order)
# To Do: sub genome_info has redundant code to extract codon => subroutine
# Perl modules for Graphs: https://metacpan.org/release/Graph
#

use strict;
use Bio::AlignIO;
use Data::Dumper;
use Graph;
use Graph::Writer::Dot;
use Graph::Traversal::DFS;
use Bio::SeqFeature::Generic;
use JSON;
use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO;
use Bio::Tools::CodonTable;
my $myCodonTable   = Bio::Tools::CodonTable->new();
#my $writer = Graph::Writer::Dot->new();
#die "Usage: $0 --hap <fasta> --vcf <vcf> --genome <fasta>\n" unless;
# Usage: perl hapnet.pl --genome ref.gb --vcf snps2-ref-03-19.vcf --hap imputed.aln --impute-log impute.log > net.json
my %options;
my $outEPI = 'EPI_ISL_402131';
my $rootEPI = 'EPI_ISL_406801'; # determined after replicated runs (H129, H2, H4, H91, ~25% each)
my $orfShift = 13468; # orf1ab reading frame shift
#my $impute_aln = 'imputed.aln';
#my $impute_log = 'impute.log';
#my $gb = 'ref.gb';
GetOptions (
    \%options,
    "hap=s", # required
    "vcf=s", # required
    "genome=s", # required
    "output=s",
    "impute-log=s", # required
    "help",
    "debug|d",
    "majority-parent=s", # read pa-cts.long file
    "add-root", # force root edge
#    "no-uniq-seq", # do not consolidate into uniq STs (for bootstrap analysis, with STs as input)
    ) or pod2usage(2);

pod2usage(1) if $options{'help'};
#print Dumper(\%options); exit;
#####################################
# Read Files
#######################################
my %hapST;
my %STs;

# Read impute log file
open IMP, "<", $options{'impute-log'} || die "need impute.log file\n";
my @sites_bound;
my $outST; 
my $rootST; 
while(<IMP>) {
    chmod;
    if (/^(\d+)\s+(\d+)$/) {
	@sites_bound = ($1, $2);
    }

    if (/^ST(\d+)\s+(EPI_ISL_\d+).*/) {
	my ($st, $hap) = ("H" . $1, $2); # assign Haplotype numbers
	$outST = $st if $hap eq $outEPI;
	$rootST = $st if $hap eq $rootEPI;
	if ($hapST{$st}) {
	    my $ref = $hapST{$st};
	    push @$ref, $hap;
	} else { # not seen
	    $hapST{$st} = [ ($hap) ]
	}
    }
}
close IMP;
print Dumper(\%hapST) if $options{'debug'}; 
die "root ST not found\n" unless $rootST;

#exit;

# Read VCF file
my @pos;
open VCF, "<", $options{'vcf'} || die "need a pre-imputed vcf file\n";;
my $ct_site = 0;
while(<VCF>){
    chomp;
    next unless /^NC_/;
    $ct_site ++;
    next unless $ct_site >= $sites_bound[0] && $ct_site <= $sites_bound[1]; # remove ends (based on impute.log file)
    my @a = split;
    push @pos, $a[1];
}
close VCF;
#print scalar @pos, "\n"; exit;

my $gIN = Bio::SeqIO->new(-file => $options{'genome'}, -format => 'genbank');
my $ge = $gIN->next_seq();

# Read CSV data & create genome features
my @gene_feats;
while (<DATA>) {
    chomp;
    next unless /^NC_/;
    my ($chr, $start, $end, $gid, $strand, $rest) = split; 
    push @gene_feats, Bio::SeqFeature::Generic->new(-seq_id => $gid, -start => $start, -end => $end, -strand => $strand eq 'plus' ? '+':'-');
}

print Dumper(\@gene_feats) if $options{'debug'}; 
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
print Dumper(\@igs_feats) if $options{'debug'};
my @utr_feats;
push @utr_feats, Bio::SeqFeature::Generic->new(-seq_id => '5-UTR', 
					       -start => 1, 
					       -end => $gene_feats[0]->start -1, 
					       -strand => 1);

push @utr_feats, Bio::SeqFeature::Generic->new(-seq_id => '3-UTR', 
					       -start => $gene_feats[-1]->end + 1, 
					       -end => 29903, # end of NC_045512.2 
					       -strand => 1);


print Dumper(\@utr_feats) if $options{'debug'};

#############################
# Make a graph, calculate MST & polorize with root
#############################
my %edge_diff;
my $alnLen;
# Read alignments
my $in = Bio::AlignIO->new(-file=>$options{'hap'}, -format=>'clustalw');
my $ct = 0; # could be multiple alignments (e.g., bootstrap)
while (my $aln = $in->next_aln()) {
    $ct++;
    $alnLen = $aln->length();
    foreach my $seq ($aln->each_seq) {
	my $id = $seq->display_id();
	$id =~ s/ST/H/; # hap ids, H1, H2, etc
	$seq->display_id($id);
	$STs{$id} = $seq; # hap seq obj
    }

    my (@nodelist, @edgelist);
    my $graph = &make_complete_graph();
    my $mstg = $graph->MST_Kruskal(); # Krustal's MST Algorithm

    if ($options{'add-root'}) {
    # the following enforce root to $rootST;
	my @outE = $mstg->edges_from($outST); # remove all edges from & to outST
	foreach my $e (@outE) {
#	    my ($u, $v) = ($e->[0], $e->[1]);
#	    my $other = $u eq $outST ? $v : $u;
#	    my @vs = $
	    $mstg->delete_edge($e);
	}
	$mstg->add_edge($outST, $rootST); # create a new edge to attach root with outgroup	
    }
    
    &attach_attributes_and_polarize($mstg, \@nodelist, \@edgelist);
    print Dumper($mstg) if $options{'debug'};

#############################
# Write outputs
##############################

    open JS, ">", "net-$ct.json";
    print JS to_json([\@nodelist, \@edgelist]); 
    close JS;
    &print_node_list(\@nodelist, $ct); 
    &print_edge_list($mstg, $ct); 
}
exit;

#####################################
# subroutines
###################################

sub make_complete_graph {
    my $g = Graph->new(undirected => 1);
    my @sts = sort keys %STs;
    for (my $i=0; $i<$#sts; $i++) { # build a complete graph from edges
	my $idI = $STs{$sts[$i]}->display_id();
	for (my $j=$i+1; $j<=$#sts; $j++) {
	    my $idJ = $STs{$sts[$j]}->display_id();
	    $g->add_edge($idI, $idJ);
	    my ($ref1, $ref2) = &comp_seq($STs{$sts[$i]}, $STs{$sts[$j]}); 
#	print Dumper(\@diffs);
	    $edge_diff{$idI}{$idJ} = $ref1; # save for later use
	    $edge_diff{$idJ}{$idI} = $ref2;
	    $g->set_edge_weight($idI, $idJ, scalar @$ref1); # seq diff as weight (to be minimized in MST)
	}
    }
    return $g;
}

sub attach_attributes_and_polarize {
    my $g = shift;
    my $refNode = shift;
    my $refEdge = shift;
# 1st traversal
    my $attach_haps = sub {
	my ($v, $self) = @_;
	my $ref = $hapST{$v};
	$g->set_vertex_attribute($v, 'haps', $hapST{$v});
	$g->set_vertex_attribute($v, 'id', $v);
	push @$refNode, $g->get_vertex_attributes($v);
    };

    my $traversal = Graph::Traversal::DFS->new(
	$g,
	pre => $attach_haps # process each vertex by pre-order traversal
	);
    $traversal->dfs();
    print Dumper($refNode) if $options{'debug'};

# 2nd traversal
# Call back function to process each edge using an anonymous sub
    my $polarize_an_edge = sub {
	my ($u, $v, $self) = @_;
	if ($u eq $outST) { warn "root position:", $u, "=>", $v, "\n"}
	$g->set_edge_attribute($u, $v, 'change', $edge_diff{$u}{$v});
	$g->set_edge_attribute($u, $v, 'source', $u);
	$g->set_edge_attribute($u, $v, 'target', $v);
	$g->set_edge_attribute($u, $v, 'id', join "-", ($u, $v));
	push @$refEdge, $g->get_edge_attributes($u, $v);
    };

    $traversal = Graph::Traversal::DFS->new(
	$g,
	start => $outST,
	pre_edge => $polarize_an_edge
	);
    $traversal->dfs();
}

sub print_node_list {
    my $refNode = shift;
    my $i = shift;
    open ND, ">", "nodes-$i.tsv";
    foreach my $nd (@$refNode) {
	my $ref  = $nd->{haps};
	foreach my $iso (@$ref) {
	    print ND $nd->{id}, "\t", $iso, "\n";
	}	
    }
    close ND;
}


sub print_edge_list {
    my $g = shift;
    my $i = shift;
    open EG, ">", "edges-$i.tsv";
    foreach my $e ($g->edges) {
	my ($u, $v) = ($e->[0], $e->[1]);
	my $from = $g->get_edge_attribute($u, $v, 'source');
	my $to = $g->get_edge_attribute($u, $v, 'target');	
	my $refFrom = $hapST{$u}; # print Dumper($refFrom);
	my @a = @$refFrom; @a = sort @a;
	my $hapFrom = shift @a; # take the first EPI as node ID 
	my $refTo = $hapST{$v}; # print Dumper($refTo);
	my @b = @$refTo; @b = sort @b;
	my $hapTo = shift @b; 
	my $ref = $g->get_edge_attribute($u, $v, 'change');
	foreach my $mut (@$ref) {
	    my $syn_nonsyn = 'NA';
	    if ($mut->{src_aa}) {
		$syn_nonsyn = $mut->{src_aa} eq $mut->{tgt_aa} ? 1 : 0;
	    }
	    print EG join "\t", ($from, $to, 
			      $mut->{site}, $mut->{label}, 
			      $mut->{src_base}, $mut->{tgt_base}, 
			      $mut->{src_codon} || 'NA', $mut->{tgt_codon} || 'NA',
			      $mut->{src_aa} || 'NA', $mut->{tgt_aa} || 'NA',
			      $syn_nonsyn
	    );
	    print EG "\n";
	}
    }
    close EG;
}

sub comp_seq { # count only definitive changes (assuming no change at positions with missing values)
    my ($i, $j) = @_;
    my @str1 = split '', $i->seq();
    my @str2 = split '', $j->seq();
    my $id1 = $i->display_id();
    my $id2 = $j->display_id();
    my (@diff_pos, @diff_rev);
    for (my $k=0; $k<$alnLen; $k++) {
	next unless $str1[$k] =~ /[atcg]/i; # i is missing, skip
	next unless $str2[$k] =~ /[atcg]/i; # j is missing, skip
	next if $str1[$k] eq $str2[$k]; # no diff, skip
	my $ref = &genome_info($pos[$k], $i->subseq($k+1,$k+1), $j->subseq($k+1,$k+1));
	# save two items for polorize changes
	push @diff_pos, {'site' => $pos[$k], 
			    'src_base' => $i->subseq($k+1,$k+1), 
			    'tgt_base' => $j->subseq($k+1,$k+1),
			    'label' => $ref->[0],
			    'pos' => $ref->[2],
			    'cd_pos' => $ref->[3] || undef,
			    'src_codon' => lc($ref->[4]) || undef,
			    'tgt_codon' => lc($ref->[5]) || undef,
			    'src_aa' => $ref->[6] || undef,
			    'tgt_aa' => $ref->[7] || undef,
			    'src_id' => $id1,
			    'tgt_id' => $id2,
	};
	
	push @diff_rev, {'site' => $pos[$k], 
			     'src_base' => $j->subseq($k+1,$k+1), 
			     'tgt_base' => $i->subseq($k+1,$k+1),
			     'label' => $ref->[0],
			     'pos' => $ref->[2],
			     'cd_pos' => $ref->[3] || undef,
			     'tgt_codon' => lc($ref->[4]) || undef,
			     'src_codon' => lc($ref->[5]) || undef,
			     'tgt_aa' => $ref->[6] || undef,
			     'src_aa' => $ref->[7] || undef,
			     'tgt_id' => $id1,
			     'src_id' => $id2,
	} ;
    } 
    return (\@diff_pos, \@diff_rev);
}

sub genome_info {
    my ($site, $base1, $base2) = @_;
    my $feat = &genome_region($site);
    if ($feat->seq_id =~ /UTR/ || $feat->seq_id =~ /igs/) { # IGS or UTR site
	return [ $feat->seq_id,  $feat->length, $site - $feat->start + 1 ];
    } elsif ($feat->seq_id eq 'orf1ab') { # ORF shift
	my $dna = Bio::Seq->new(-id=>'frag', -seq=>$ge->subseq($feat->start, $site));
	my $dnaLen = $dna->length();
	my $mod;
	my $codonPos;
	my ($codon1, $codon2, $aa1, $aa2);
	if ($site < $orfShift) { # before the shift site
	    $mod = $dnaLen % 3;
	    if ($mod == 1) { # 1st codon position
		$codon1 = $base1 . $ge->subseq($site+1, $site+2); 
		$codon2 = $base2 . $ge->subseq($site+1, $site+2);
		$codonPos = 1;
	    } elsif ($mod == 2) {
		$codon1 = $ge->subseq($site-1, $site-1) . $base1 . $ge->subseq($site+1, $site+1); 
		$codon2 = $ge->subseq($site-1, $site-1) . $base2 . $ge->subseq($site+1, $site+1);
		$codonPos = 2; 
	    } else { # 3rd codon position
		$codon1 = $ge->subseq($site-2, $site-1) . $base1; 
		$codon2 = $ge->subseq($site-2, $site-1) . $base2; 
		$codonPos = 3;
	    }
	    $aa1 = $myCodonTable->translate($codon1);
	    $aa2 = $myCodonTable->translate($codon2);
	    return [ $feat->seq_id,  $feat->length, $site - $feat->start + 1, $codonPos, $codon1, $codon2, $aa1, $aa2 ];
	} elsif ($site == $orfShift) {
	    $mod = 0; # at the orf shift site
	    $codon1 = $ge->subseq($site-2, $site-1) . $base1; 
	    $codon2 = $ge->subseq($site-2, $site-1) . $base2; 
	    $codonPos = 3;
	    $aa1 = $myCodonTable->translate($codon1);
	    $aa2 = $myCodonTable->translate($codon2);
	    return [ $feat->seq_id,  $feat->length, $site - $feat->start + 1, $codonPos, $codon1, $codon2, $aa1, $aa2 ];
	} else { # after the shift site
	    $mod = ($dnaLen + 1) % 3; # 1->2; 2->3; 3->1
	    if ($mod == 1) { # 1st codon position
		$codon1 = $base1 . $ge->subseq($site+1, $site+2); 
		$codon2 = $base2 . $ge->subseq($site+1, $site+2);
		$codonPos = 1;
	    } elsif ($mod == 2) {
		$codon1 = $ge->subseq($site-1, $site-1) . $base1 . $ge->subseq($site+1, $site+1); 
		$codon2 = $ge->subseq($site-1, $site-1) . $base2 . $ge->subseq($site+1, $site+1);
		$codonPos = 2; 
	    } else { # 3rd codon position
		$codon1 = $ge->subseq($site-2, $site-1) . $base1; 
		$codon2 = $ge->subseq($site-2, $site-1) . $base2; 
		$codonPos = 3;
	    }
	    $aa1 = $myCodonTable->translate($codon1);
	    $aa2 = $myCodonTable->translate($codon2);
	    return [ $feat->seq_id,  $feat->length, $site - $feat->start + 1, $codonPos, $codon1, $codon2, $aa1, $aa2 ];
	}
    } else { # other ORFs
	my $dna = Bio::Seq->new(-id=>'frag', -seq=>$ge->subseq($feat->start, $site));
	my $dnaLen = $dna->length();
	my $mod = $dnaLen % 3;
	my $codonPos;
	my ($codon1, $codon2, $aa1, $aa2);
	if ($mod == 1) { # 1st codon position
	    $codon1 = $base1 . $ge->subseq($site+1, $site+2); 
	    $codon2 = $base2 . $ge->subseq($site+1, $site+2);
	    $codonPos = 1;
	} elsif ($mod == 2) {
	    $codon1 = $ge->subseq($site-1, $site-1) . $base1 . $ge->subseq($site+1, $site+1); 
	    $codon2 = $ge->subseq($site-1, $site-1) . $base2 . $ge->subseq($site+1, $site+1);
	    $codonPos = 2; 
	} else { # 3rd codon position
	    $codon1 = $ge->subseq($site-2, $site-1) . $base1; 
	    $codon2 = $ge->subseq($site-2, $site-1) . $base2; 
	    $codonPos = 3;
	}
	$aa1 = $myCodonTable->translate($codon1);
	$aa2 = $myCodonTable->translate($codon2);
	return [ $feat->seq_id,  $feat->length, $site - $feat->start + 1, $codonPos, $codon1, $codon2, $aa1, $aa2 ];
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
    die "feat not found: $site\n" unless $feat->seq_id;
#    print $feat->seq_id, ":", $feat->start, ":", $feat->end, "\n"; 
    return $feat;
}


__DATA__
Accession       Start   Stop    Gene symbol     Strand  NCBI Gene ID    Name
NC_045512.2     266     21555   orf1ab  plus    43740578
NC_045512.2     21563   25384   S       plus    43740568
NC_045512.2     25393   26220   ORF3a   plus    43740569
NC_045512.2     26245   26472   E       plus    43740570
NC_045512.2     26523   27191   M       plus    43740571
NC_045512.2     27202   27387   ORF6    plus    43740572
NC_045512.2     27394   27759   ORF7a   plus    43740573
NC_045512.2     27756   27887   ORF7b   plus    43740574
NC_045512.2     27894   28259   ORF8    plus    43740577
NC_045512.2     28274   29533   N       plus    43740575
NC_045512.2     29558   29674   ORF10   plus    43740576

