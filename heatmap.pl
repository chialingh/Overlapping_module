#!/usr/bin/perl

use strict;
use warnings;

my %NeXO;
open(IN, "/home/clhuang/lab/databases/NeXO/NeXO_terms.txt");
while(<IN>){
	chomp $_;
	my ($name, $genes) = split(/\s=\s/, $_);
	my @gg = split(/\+/, $genes);
	foreach my $g(@gg){
		$g =~ s/\s//g;
		$NeXO{$name}{$g} = 1;
	}
}
close IN;

my %Modules;
open(IN, "$ARGV[0]");
while(<IN>){
	chomp $_;
	my @line = split(/\t/, $_);
	my $name = shift @line;
	foreach my $g(@line){
		$g =~ s/\s//g;
		$Modules{$name}{$g} = 1;
	}
}
close IN;

my %aa;
foreach my $m1(keys %Modules){
	foreach my $m2(keys %NeXO){
		my %intersection;
		foreach(keys %{$Modules{$m1}}){
			$intersection{$_} = ${$Modules{$m1}}{$_} if exists ${$NeXO{$m2}}{$_};
		}
		my %union;
		@union{keys %{$Modules{$m1}}, keys %{$NeXO{$m2}}} = ();
		my $ins = scalar keys %intersection;
		my $uni = scalar keys %union;
		my $v = $ins*2/($uni+$ins);
		if($v > 0.3){
			$aa{$m1} = 1;
		}
	}
}

my $i = scalar keys %aa;
print "$i\n";
