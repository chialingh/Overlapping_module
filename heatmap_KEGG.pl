#!/usr/bin/perl

use strict;
use warnings;

# read modules

my @modules;
open(IN, "results/Merged_union_0p2_0.1.txt");
while(<IN>){
	chomp $_;
	my @line = split(/\t/, $_);
	my $name = shift @line;
	my $rec = {};
	foreach my $g(@line){
		$rec->{$g} = ();
	}
	push(@modules, $rec);
}
close IN;

# read KEGG pathways
open(IN, "/home/clhuang/lab/databases/KEGG/PathGene/path1.list");
while(my $path=<IN>){
	chomp $path;
	my %KEGG;
	open(INP, "/home/clhuang/lab/databases/KEGG/PathGene/$path");
	while(<INP>){
		chomp $_;
		$KEGG{$_} = ();
	}
	close INP;

	# compare modules and KEGG pathway
	foreach my $href(@modules){
		my %intersection;
		foreach(keys %KEGG){
			$intersection{$_} = $KEGG{$_} if exists ${$href}{$_};
		}
		my $sn = scalar keys %intersection;
		my $n1 = scalar keys %$href;
		my $n2 = scalar keys %KEGG;
		my $p = $sn/($n1+$n2-$sn);
		print "$sn\t";
	}
	print "\n";
}
close IN;
