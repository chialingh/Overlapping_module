#!/usr/bin/perl

use strict;
use warnings;


my $i = 0;
my $j = 0;
my @modules;
open(IN, "aa2");
while(<IN>){
	chomp $_;
	if($i%2){
		my $rec = {};
		my @line = split(/\t/, $_);
		foreach my $g(@line){
			$rec->{$g} = ();
		}
		push(@modules, $rec);
	}
	$i = $i + 1;
	$j = $j + 1;
}
close IN;

foreach my $href1(@modules){
	foreach my $href2(@modules){
		my %intersection;
		foreach(keys %$href1){
			$intersection{$_} = ${$href1}{$_} if exists ${$href2}{$_};
		}
		my $sn = scalar keys %intersection;
		my $n1 = scalar keys %$href1;
		my $n2 = scalar keys %$href2;
		my $p = $sn/($n1+$n2-$sn);
		print "$sn\t";
	}
	print "\n";
}
