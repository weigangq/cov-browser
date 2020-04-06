#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

my %seen;
while(<>) {
    chomp;
    next unless /^\s+(\d+)\s+(H\d+)-(H\d+)$/;
    my ($ct, $parent, $child) = ($1, $2, $3);
    $seen{$child}->{$parent} = $ct;
}

#print Dumper(\%seen);

foreach my $child (keys %seen) {
    my %ref = %{$seen{$child}};
    my @pas = sort { $ref{$b} <=> $ref{$a}} keys %ref;
    foreach my $pa (@pas) {
	print join "\t", ($child, $pa, $seen{$child}->{$pa});
	print "\n";
    }
#    print $child;
#    foreach (@pas) {
#	print "\t", ($_, "|", $seen{$child}->{$_});
#    }
#    print "\n";
}
exit;
