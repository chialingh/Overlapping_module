#!/usr/bin/perl
# take mutiple seeds
# Do not check if input seeds are module
# simulated annealing
# FLN

use strict;

my $file1 = $ARGV[0]; # seeds
my $theda = $ARGV[1];
#my $file2 = "FLN.data.human.txt"; # FLN
#my $NN = 21659;
my $file2 = "random_small_FLN.txt";
my $NN = 11015;

my ($f, $apx) = split(/\./, $file1);

my @seeds;
my %FLN;
my %module;

# read FLN
open(IN, "/home/clhuang/lab/FLN/$file2");
while(<IN>){
	chomp $_;
	my @line = split(/\t/, $_);
	$FLN{$line[0]}{$line[1]} = abs($line[2]);
	$FLN{$line[1]}{$line[0]} = abs($line[2]);
}
close IN;

my %module; # current module
my %nodes; # looked nodes
my @all_ns;

# read seeds
#open(IN, "/home/clhuang/lab/causal_modules/results/2012NOV06/seeds/$file1");
open(IN, "seeds1.txt");
while(<IN>){
	chomp $_;
	$nodes{$_} = (); # looked nodes
	$module{$_} = (); # first module
}
close IN;

my $e = 2.71;
my $lamda = 0.9;
my $Sd = cal_mod_score(\%module);

for(my $T = 70; $T > 0.00005; $T = $lamda*$T){
	my %new_ns; # neighbors of first layer nodes
	foreach my $g1(keys %module){
		foreach my $g2(keys %{$FLN{$g1}}){
			$new_ns{$g2} = ();
		}
	}
	my @new_ns_array  = keys %new_ns;
	my $g2_size = scalar @new_ns_array;
	my $i = int(rand($g2_size));
	my $new_node = $new_ns_array[$i];
	my %new_module = %module;
	$new_module{$new_node} = ();

	my $Sd_new = cal_mod_score(\%new_module);
	if($Sd_new <= $Sd){
		$Sd = $Sd_new;
		%module = %new_module;
	}else{
		if(rand(1) <= $e**(-1*($Sd_new-$Sd)/$T)){
			$Sd = $Sd_new;
			%module = %new_module;
		}
	}
}

# print results
foreach(keys %module){
	print "$_\n";
}

sub cal_mod_score{

	my $mod_ref = shift;
	my %module = %$mod_ref; # module

	my $Si1 = 0;
	my $So1 = 0;
	my $e1 = 0;
	my $e2 = 0;
	foreach my $g1(keys %module){ # for every node in current module
		foreach my $g2(keys %{$FLN{$g1}}){ # look for their neighbors on FLN
			if(exists $module{$g2}){ # if a neighbor is inside module, calculate Si
				$Si1 = $Si1 + $FLN{$g1}{$g2};
				$e1 = $e1 + 1;
			}else{ # if a neighbor is outside module, calculate So
				$So1 = $So1 + $FLN{$g1}{$g2};
				$e2 = $e2 + 1;
			}
		}
	}


	# in_score1 and out_score1 are scores of current module
	my $in_score1 = 0;
	my $out_score1 = 0;
	$Si1 = $Si1/2;
	$e1 = $e1/2;
	if($e1 > 0){
		$in_score1 = $Si1/$e1;
	}
	if($e2 > 0){
		$out_score1 = $So1/$e2;
	}

	my $Sd1 = (-1)*($in_score1 - $out_score1); # score of current module

	return($Sd1);
}

