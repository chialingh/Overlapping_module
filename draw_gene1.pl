#!/usr/bin/perl

use strict;
#use warnings;

my @FLN_gene;
open(IN, "/home/clhuang/lab/FLN/FLN.sce.gene");
#open(IN, "/home/clhuang/lab/databases/NeXO/nbt.ent");
while(<IN>){
        chomp $_;
#	my ($a, $b) = split(/\t/, $_);
        push(@FLN_gene, $_);
}
close IN;

my $N = scalar @FLN_gene;

my %FLN;
open(IN, "/home/clhuang/lab/FLN/FLN.data.sce.txt");
#open(IN, "/home/clhuang/lab/databases/NeXO/aa");
while(<IN>){
	chomp $_;
	my ($g1, $g2, $w) = split(/\t/, $_);
#	my ($g1, $g2) = split(/\t/, $_);
	$FLN{$g1}{$g2} = 1;
	$FLN{$g2}{$g1} = 1;
}
close IN;

my $i = 1;
while($i <= 10000){
	my @seeds;
	my %seeds_hash;
	my %pool;
	$seeds[0] = $FLN_gene[int(rand($N))];
	$seeds_hash{$seeds[0]} = 1;
	for(my $j=0; $j<=3; $j++){
		foreach(keys %{$FLN{$seeds[$j]}}){
			$pool{$_} = ();
		}
		my @pool_array = keys %pool;
		last if scalar @pool_array < 1;
		my $np = scalar @pool_array;
		my $d = int(rand($np));
		push(@seeds, $pool_array[$d]) unless $seeds_hash{$pool_array[$d]};
		$seeds_hash{$pool_array[$d]} = 1;
	}

	my $z = scalar @seeds;

	if($z == 5){
		open(OUT, ">SCE/seeds1/$i.sds");
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

