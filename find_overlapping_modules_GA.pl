#!/usr/bin/perl
# take mutiple seeds
# Do not check if input seeds are module
# simulated annealing
# FLN

use strict;
use warnings;

my $file1 = $ARGV[0]; # seeds
#my $theda = $ARGV[1];
#my $file2 = "FLN.data.human.txt"; # FLN
#my $NN = 21659;
#my $file2 = "test_FLN.txt";
my $file2 = "FLN.data.sce.txt";
my $NN = 5469;

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

my %nodes; # looked nodes
my @all_ns;

# read seeds
#open(IN, "/home/clhuang/lab/causal_modules/results/2012NOV06/seeds/$file1");
open(IN, "$file1");
while(<IN>){
	chomp $_;
	$nodes{$_} = (); # looked nodes
	$module{$_} = (); # first module
}
close IN;

my $e = 2.71;
my $lamda = 0.99;
my $Sd = cal_mod_score(\%module);

for(my $T = 100; $T > 0.0005; $T = $lamda*$T){
	my %new_ns; # neighbors of first layer nodes
	foreach my $g1(keys %module){
		foreach my $g2(keys %{$FLN{$g1}}){
#			my %temp_module = %module;
#			$temp_module{$g2} = ();
#			my $temp_Sd = cal_mod_score(\%temp_module);
#			if(($temp_Sd - $Sd) <= $theda*(-1)){
				$new_ns{$g2} = ();
#			}
		}
	}
	my $j = 0;
	for(my $i = 1; $i <= 200; $i++){
		my $aa = scalar keys %module;
		exit if $aa == 0;
		# randomly pick one node from neighbors and add it into the module
		my @new_ns_array  = keys %new_ns;
		my $g2_size = scalar @new_ns_array;
		my $j = int(rand($g2_size));
		my $new_node = $new_ns_array[$j];
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
		# randomly pick one node from the module and toss it
		my @module_nodes_array = keys %module;
		my $mod_size = scalar @module_nodes_array;
		$j = int(rand($mod_size));
		my $toss_node = $module_nodes_array[$j];
		delete $new_module{$toss_node};
		$Sd_new = cal_mod_score(\%new_module);
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
}

my $aa = scalar keys %module;
# print results
print "$aa\t$Sd\n";
my @ee = keys %module;
my $bb = shift @ee;
print "$bb";
foreach(@ee){
	print "\t$_";
}
print "\n";

sub cal_mod_score{

	my $mod_ref = shift;
	# my %module = %$mod_ref; # module

	my $Si1 = 0;
	my $So1 = 0;
	my $e1 = 0;
	my $e2 = 0;
	foreach my $g1(keys %$mod_ref){ # for every node in current module
		foreach my $g2(keys %{$FLN{$g1}}){ # look for their neighbors on FLN
			if(exists ${$mod_ref}{$g2}){ # if a neighbor is inside module, calculate Si
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

#	my $Sd1 = $in_score1 - $out_score1; # score of current module
	my $Sd1 = $out_score1 - $in_score1; # score of current module

	return($Sd1);
}

