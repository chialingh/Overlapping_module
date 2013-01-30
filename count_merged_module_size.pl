#!/usr/bin/perl

use strict;

open(IN, "$ARGV[0]");
while(<IN>){
	chomp $_;
	my @line = split(/\t/, $_);
	my $name = shift @line;
	my $size = scalar @line;
	print "$name\t$size\n";
}
close IN;
