#!/usr/bin/perl

use strict;

my @FLN_gene;
#open(IN, "/home/clhuang/lab/FLN/FLN.human.gene");
open(IN, "/home/clhuang/lab/databases/NeXO/nbt.ent");
while(<IN>){
        chomp $_;
	my ($a, $b) = split(/\t/, $_);
        push(@FLN_gene, $b);
}
close IN;

my %FLN;
#open(IN, "/home/clhuang/lab/FLN/FLN.data.human.txt");
open(IN, "/home/clhuang/lab/databases/NeXO/aa");
while(<IN>){
	chomp $_;
#	my ($g1, $g2, $w) = split(/\t/, $_);
	my ($g1, $g2) = split(/\t/, $_);
	$FLN{$g1}{$g2} = 1;
	$FLN{$g2}{$g1} = 1;
}
close IN;

my $i = 1;
while($i <= 10000){
	my @rand_seeds;
	my @seeds;
	my $rnd_ref = fisher_yates_shuffle( \@FLN_gene);
	my @rand_seeds = @$rnd_ref;
	my @seeds = @rand_seeds[0..4];

	my %sds;

	# use direct connected genes as seeds
	foreach my $s1(@seeds){
		foreach my $s2(@seeds){
			if($FLN{$s1}{$s2}){
				$sds{$s1} = 1;
			}
		}
	}

	my $z = scalar keys %sds;

	if($z == 5){
		open(OUT, ">NeXO/seeds/$i.sds");
		print OUT "$_\n" foreach @seeds;
		close OUT;
		$i = $i + 1;
	}
}

sub fisher_yates_shuffle{
    my $array = shift;
    my $i;
    for ($i = @$array; --$i; ) {
        my $j = int rand ($i+1);
        next if $i == $j;
        @$array[$i,$j] = @$array[$j,$i];
    }
    my @array = @$array;
    return \@array;
}

