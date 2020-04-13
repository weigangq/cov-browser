#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;
use JSON;
use Bio::TreeIO;
use Getopt::Long;

my %opts;
GetOptions (
    \%opts,
    "ladder",
    "debug|d",
    );
die "Usage: $0 <network dnd in NEWICK>\n" unless @ARGV == 1;
my $filename = shift;

#read tree:
my $in = Bio::TreeIO->new(-file=>$filename, -format=>'newick');
my $tree = $in->next_tree();
#print Dumper($tree); exit;
#coord of nodes:
my $root = $tree->get_root_node();
my @nodes = $tree->get_nodes();
my %desOTUs; 
my %desAll;# counts of all desc nodes for each internal node
foreach my $nd (@nodes) { 
    my $bl = $nd->branch_length() || 0; 
    $nd->branch_length($bl);
    my $decCts = 0;
    my @decLeaf;
    if ($nd->is_Leaf) {
	$desAll{$nd} = $decCts;
	next;
    }

    foreach ($nd->get_all_Descendents) {
	$decCts++;
	push @decLeaf, $_->id() if $_->is_Leaf;
    }
    $desAll{$nd} = $decCts;
    $desOTUs{$nd} = \@decLeaf;
}

my $sp_l = 1; # space between OTUs
my $ymin = 0;
my (%idY, @nodes_js, @leafs);

if ($opts{'ladder'}) {
    &ladder_ycoord_ref($root, \$ymin); # assign y coord for all nodes recursively
} else {
    &assign_leaf_ycoord_ref($root, \$ymin); # assign y coord for leaf nodes (used to start from refnode; now start from root)
    &assign_inode_ycoord($root); # assign y coord for internal nodes (few)
}

foreach my $node (@nodes) {
    $node->{xcoord} = &distance_to_root($node);
    my %n;
    $n{xcoord} = $node->{xcoord};
    $n{ycoord} = $node->{ycoord};
    $n{branch_length} = $node->branch_length || 0;
    $n{id} = $node->id();
    $idY{$node->id()} = $node->{ycoord};

    if ($node->is_Leaf()) {
	$n{is_Leaf} = 1;
    } else {
	my @tmp;
	push @tmp, $_->{'ycoord'} foreach $node->each_Descendent();
	my @sorted = sort { $a <=> $b } @tmp;
	$n{descendent_ycoord} = [$sorted[0], $sorted[-1]];
#	$n{descendent_ycoord} = \@sorted;
	
#	@leafs = ();
#	&get_leafs($node);
#	my @ll = @leafs;
	$n{leaves} = $desOTUs{$node};
	$n{parent_is_root} = 1 if $node ne $root && $node->ancestor() eq $root;
    }
    push @nodes_js, \%n
}

my $w_tree = ($root->branch_length || 0) + $root->height;

my $unit=sprintf "%.0e", $w_tree/10;
$unit *= 1;
my $ll = substr($unit, -1);
$unit *= $ll==3 ? 2/3 : ($ll==4 ? 5/4 : ($ll==6 ? 5/6 : ($ll==7 ? 5/7 : ($ll==8 ? 10/8 : ($ll==9 ? 10/9 : 1)))));

#output @ids:
my @ids = sort {$idY{$a} <=> $idY{$b}} keys %idY;
#print STDERR join ",", @ids;
#print STDERR "\n";

#output tree:
my $output;
$output->{ids} = \@ids;
$output->{nodes} = \@nodes_js;
$output->{w_tree} = $w_tree;
$output->{unit} = $unit;
#print Dumper(\@nodes_js);
print to_json($output);

exit;

sub ladder_ycoord_ref {
    my ($node, $yposref) = @_;
    $node->{'ycoord'} = $$yposref; 
    return if $node->is_Leaf();
    my @desc = sort {$desAll{$a} <=> $desAll{$b}} ( $node->each_Descendent() ); # sort from less to more descendant nodes

    if (@desc == 1) { # single child & an internal node
	$desc[0]->{'ycoord'} = $$yposref; # don't increment
	&ladder_ycoord_ref ($desc[0], $yposref); # recurse on internal-node child
    } else { # 2 or more children
	for (my $i=0; $i<=$#desc; $i++ ) {
	    $desc[$i]->{'ycoord'} = $$yposref; # first child at the same y ... 
	    &ladder_ycoord_ref ($desc[$i], $yposref); # recurse on internal-node child
	    $$yposref += $sp_l unless $i == $#desc  # increment only when more than one child unless it is the last
	}
    }
}


sub assign_leaf_ycoord_ref {
    my ($node, $yposref) = @_;
    if ($node->is_Leaf()) { 
	$node->{'ycoord'} = $$yposref; 
	warn "Hit OTU: ", $node->id(), " at $$yposref\n";
	$$yposref += $sp_l;
	return;
    } # increment y only when OTU

    for my $child ( $node->each_Descendent() ) { 
	&assign_ycoord_ref ($child, $yposref)
    }
}

sub assign_inode_ycoord {
  my $node = shift;
  my @tmp;
  return if $node->is_Leaf;
  for my $child ( $node->each_Descendent() ) {
    assign_inode_ycoord($child) unless $child->{'ycoord'};
    push @tmp, $child->{'ycoord'};
  }
  if (@tmp > 1) { # at least two desc
      my @sorted = sort { $a <=> $b } @tmp;
      $node->{'ycoord'} = sprintf "%.4f", $sorted[0] + 1 / 2 * ( $sorted[-1] - $sorted[0] );
  } else { # only one desc
      $node->{'ycoord'} = $tmp[0];
  }
  warn "innode: ", $node->id() || $node->internal_id(), " at ", $node->{'ycoord'}, "\n";
}

sub distance_to_root {
  my $node = shift;
  my $dist = $node->branch_length() || 0;
  my $parent = $node->ancestor();
  $dist += &distance_to_root($parent) if defined $parent;
  return $dist
}



