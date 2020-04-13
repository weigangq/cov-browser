#!/usr/bin/perl -w
use warnings;
use strict;
use Data::Dumper;
use JSON;
use Bio::TreeIO;

die "Usage: $0 <treefile in NEWICK> <ref_node>\n" unless @ARGV == 2;

my $filename = shift;
my $ref = shift;

#read tree:
my $in = Bio::TreeIO->new(-file=>$filename, -format=>'newick');
my $tree = $in->next_tree();


#coord of nodes:
my $root = $tree->get_root_node();
my @nodes = $tree->get_nodes();
foreach my $nd (@nodes) { my $bl = sprintf "%.10g", $nd->branch_length() || 0; $nd->branch_length($bl) }

my $ref_node;
foreach my $node (@nodes) {
  next unless $node->is_Leaf();
  $ref_node = $node if $node->id() eq $ref;
}

my $sp_l = 1;
my $ymin = 0;
&assign_leaf_ycoord_ref($ref_node, \$ymin);

&assign_inode_ycoord($root);

my (%leaves, @nodes_js, @leafs);
foreach my $node (@nodes) {
  $node->{xcoord} = &distance_to_root($node);

  my %n;
  $n{xcoord} = $node->{xcoord};
  $n{ycoord} = $node->{ycoord};
  $n{branch_length} = $node->branch_length if $node->branch_length;

  if ($node->is_Leaf()) {
    $n{is_Leaf} = 1;
    $leaves{$node->id()} = $node->{ycoord};
    $n{id} = $node->id()
  } else {
    my @tmp;
    push @tmp, $_->{'ycoord'} foreach $node->each_Descendent();
    my @sorted = sort { $a <=> $b } @tmp;
    $n{descendent_ycoord} = [$sorted[0], $sorted[-1]];

    @leafs = ();
    get_leafs($node);
    my @ll = @leafs;
    $n{leaves} = \@ll;
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
my @ids = sort {$leaves{$a} <=> $leaves{$b}} keys %leaves;
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


sub assign_leaf_ycoord_ref {
    my ($node, $yposref) = @_;

    if ($node->is_Leaf()) { $node->{'ycoord'} = $$yposref; $$yposref += $sp_l }

    my $sis = &get_sister($node);
    &walk_otu($_, $yposref) foreach @$sis;

    my $parent = $node->ancestor();
    return if $parent eq $root;
    &assign_leaf_ycoord_ref ($parent, $yposref)
}

sub get_sister {
  my $node = shift;
  my $parent = $node->ancestor();
  my @nds;
  for my $child ( $parent->each_Descendent() ) { push @nds, $child unless $child eq $node }
  return \@nds
}

sub walk_otu {
    my ( $node, $yposref) = @_;
    if ($node->is_Leaf()) {
      $node->{'ycoord'} = $$yposref;
      $$yposref += $sp_l;
      return
    } else {
      &walk_otu($_, $yposref) for $node->each_Descendent()
    }
}

sub assign_inode_ycoord {
  my $node = shift;
  my @tmp;
  for my $child ( $node->each_Descendent() ) {
    assign_inode_ycoord($child) unless defined $child->{'ycoord'};
    push @tmp, $child->{'ycoord'}
  }
  my @sorted = sort { $a <=> $b } @tmp;
  $node->{'ycoord'} = $sorted[0] + 1 / 2 * ( $sorted[-1] - $sorted[0] )
}

sub distance_to_root {
  my $node = shift;
  my $dist = $node->branch_length() || 0;
  my $parent = $node->ancestor();
  $dist += &distance_to_root($parent) if defined $parent;
  return $dist
}

sub get_leafs {
  my $node = shift;
  for my $child ( $node->each_Descendent() ) {
    if ($child->is_Leaf) { push @leafs, $child->id() }
    else { &get_leafs($child) }
  }
}
