#!/usr/bin/perl
# take mutiple seeds
# Do not check if input seeds are module
# FLN

use strict;
use warnings;
use Math::BigInt;

my $file1 = $ARGV[0]; # seeds
my $theda = $ARGV[1];
my $NN = 5070;

my ($f, $apx) = split(/\./, $file1);

my @seeds;
my %FLN;
my %module; # current module

# read FLN
# open(IN, "/home/clhuang/lab/FLN/$file2");
open(IN, "/home/clhuang/lab/databases/NeXO/aa");
while(<IN>){
	chomp $_;
	my @line = split(/\t/, $_);
	$FLN{$line[0]}{$line[1]} = 1;
	$FLN{$line[1]}{$line[0]} = 1;
}
close IN;

my %nodes; # looked nodes

# read seeds
# open(IN, "/home/clhuang/lab/causal_modules/results/2012NOV06/seeds/$file1");
open(IN, "/home/clhuang/lab/causal_modules/results/Final/Overlapping_module/NeXO/seeds1/$file1");
while(<IN>){
	chomp $_;
	$nodes{$_} = 1; # looked nodes
	$module{$_} = 1; # first module
}
close IN;

my %new_ns; # neighbors of first layer nodes
foreach my $g1(keys %module){
	foreach my $g2(keys %{$FLN{$g1}}){
		unless($module{$g2}){ # if new neighbor g2 is not in current module
			$new_ns{$g2} = 1;
		}
	}
}

my @all_ns = keys %new_ns;

my $N = scalar keys %nodes;

while($N < $NN){ # 21659: total number of nodes in KeggFLN_all
	my ($new_mem, $nodes) = eval_node(\@all_ns, \%nodes, \%module);
	my %new_mem_hash = %$new_mem;
	my %nodes = %$nodes;
	last if scalar keys %new_mem_hash == 0;
	foreach my $m(keys %new_mem_hash){
		$module{$m} = 1;
	}


	@all_ns = ();
	my %new_ns; # neighbors of nodes in current module
	foreach my $g1(keys %module){
		foreach my $g2(keys %{$FLN{$g1}}){
			unless($module{$g2}){ # if new neighbor g2 is not in current module
				$new_ns{$g2} = 1;
			}
		}
	}
	last if scalar keys %new_ns == 0; # stop if there is no new neighbor
	@all_ns = keys %new_ns;		
}
			
# print result
# open(OUT, ">/home/clhuang/lab/causal_modules/results/2012NOV06/modules/0p15/$f.txt");
open(OUT, ">/home/clhuang/lab/causal_modules/results/Final/Overlapping_module/NeXO/modules/0p1/$f.txt");
foreach my $g(keys %module){
	print OUT "$g\n";
}
close OUT;



sub eval_node{
	my ($ns_ref, $nd_ref, $mod_ref) = @_;

	my @new_ns = @$ns_ref; # putative new neighboring nodes
	my %nd = %$nd_ref; # looked nodes
	my %mod = %$mod_ref; # current module

	my %members; # new members

	my $mod_n = scalar keys %mod;
	my $new_ns_n = scalar @new_ns;
	my $Si1 = 0;
	my $So1 = 0;
	my $e1 = factorial_p($mod_n)/factorial_p($mod_n-2);
	my $e2 = $new_ns_n * $mod_n;
	foreach my $g1(keys %mod){ # for every node in current module
		foreach my $g2(keys %{$FLN{$g1}}){ # look for their neighbors on FLN
			if($mod{$g2}){ # if a neighbor is inside module, calculate Si
				$Si1 = $Si1 + $FLN{$g1}{$g2};
				# $e1 = $e1 + 1;
			}else{ # if a neighbor is outside module, calculate So
				$So1 = $So1 + $FLN{$g1}{$g2};
				# $e2 = $e2 + 1;
			}
		}
		$nd{$g1} = 1;
	}

	last if $So1 == 0;

	# in_score1 and out_score1 are scores of current module
	$Si1 = $Si1/2;
	my $in_score1 = $Si1/$e1;
	my $out_score1 = $So1/$e2;
=head
	$e1 = $e1/2;
	if($e1 > 0){
		$in_score1 = $Si1/$e1;
	}
	if($e2 > 0){
		$out_score1 = $So1/$e2;
	}
=cut
	my $Sd1 = $in_score1 - $out_score1; # score of current module

	# for every new node g1, examine if adding g1 will increase the total score Sd of module
	for(my $i = 0; $i <= $#new_ns; $i++){
		my $Si2 = 0;
		my $So2 = 0;
		my $ee1 = factorial_p($mod_n + 1)/factorial_p($mod_n + 1 - 2);
		my $g1 = $new_ns[$i];
		$nd{$g1} = 1;
		foreach my $g2(keys %{$FLN{$g1}}){ # for every neighbor of g1
			# examine if neighbor g2 is in the current module
			if($mod{$g2}){ 
				$Si2 = $Si2 + $FLN{$g1}{$g2};
				# $ee1 = $ee1 + 1;
			}else{
				$So2 = $So2 + $FLN{$g1}{$g2};
				# $ee2 = $ee2 + 1;
			}
		}
		if($So2){
			my $ee2 = ($So2 + $new_ns_n - 1)*($mod_n + 1);
			my $in_score2 = ($Si1+$Si2)/$ee1;
			my $out_score2 = ($So1-$Si2+$So2)/$ee2;
			my $Sd2 = $in_score2 - $out_score2; # score of module including node g1
			if($Sd2 - $Sd1 > $theda){
				$members{$g1} = 1;
			}
		}else{
			$members{$g1} = 1;
		}

	}
	return (\%members, \%nd);
}

sub factorial_p{
	my $n = shift;
	my $fact = Math::BigInt->new($n);
	return $fact->bfac();
}
