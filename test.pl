#!/usr/bin/perl
# take mutiple seeds
# Do not check if input seeds are module
# simulated annealing
# FLN

use strict;
use warnings;

my $file1 = $ARGV[0]; # seeds
my $file2 = "FLN.data.sce.txt";
my $NN = 5469;

my ($f, $apx) = split(/\./, $file1);

my @seeds;
my %FLN;
my %FLN_genes;
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
open(IN, "$file1");
while(<IN>){
	chomp $_;
	$nodes{$_} = (); # looked nodes
	$module{$_} = (); # first module
}
close IN;

my $e = 2.71;
my $lamda = 0.95;
my ($Sd1, $Si1, $So1, $ei1, $eo1) = cal_mod_score(\%module);

my %new_ns; # neighbors of first layer nodes
foreach my $g1(keys %module){
	foreach my $g2(keys %{$FLN{$g1}}){
		$new_ns{$g2} = ();
	}
}

print "\n";
print "$Sd1\t\t$Si1\t\t$So1\t\t$ei1\t\t$eo1\n";
print "==============\n";
=head
		# randomly pick one node from neighbors and add it into the module
		my $aa = scalar keys %module;
		exit if $aa == 0;
		my @new_ns_array  = keys %new_ns;
		my $g2_size = scalar @new_ns_array;
		my $j = int(rand($g2_size));
		my $new_node = $new_ns_array[$j];
		my ($Sd2, $Si2, $So2, $ei2, $eo2) = cal_add_mod_score(\%module, $new_node, $Si1, $So1, $ei1, $eo1);

print "$Sd2\t\t$Si2\t\t$So2\t\t$ei2\t\t$eo2\n";
		my %new_module = %module;
		$new_module{$new_node} = ();

		($Sd2, $Si2, $So2, $ei2, $eo2) = cal_mod_score(\%new_module);

print "$Sd2\t\t$Si2\t\t$So2\t\t$ei2\t\t$eo2\n";
print "==============\n";
=cut
		my @module_nodes_array = keys %module;
		my $mod_size = scalar @module_nodes_array;
		my $j = int(rand($mod_size));
		my $toss_node = $module_nodes_array[$j];
		my ($Sd2, $Si2, $So2, $ei2, $eo2) = cal_toss_mod_score(\%module, $toss_node, $Si1, $So1, $ei1, $eo1);
print "$Sd2\t\t$Si2\t\t$So2\t\t$ei2\t\t$eo2\n";

		delete $module{$toss_node};
		($Sd2, $Si2, $So2, $ei2, $eo2) = cal_mod_score(\%module);

print "$Sd2\t\t$Si2\t\t$So2\t\t$ei2\t\t$eo2\n";
print "\n";


sub cal_mod_score{

	my $mod_ref = shift;

	my $Si = 0;
	my $So = 0;
	my $e1 = 0;
	my $e2 = 0;
	foreach my $g1(keys %$mod_ref){ # for every node in current module
		foreach my $g2(keys %{$FLN{$g1}}){ # look for their neighbors on FLN
			if(exists ${$mod_ref}{$g2}){ # if a neighbor is inside module, calculate Si
				$Si = $Si + $FLN{$g1}{$g2};
				$e1 = $e1 + 1;
			}else{ # if a neighbor is outside module, calculate So
				$So = $So + $FLN{$g1}{$g2};
				$e2 = $e2 + 1;
			}
		}
	}


	# in_score1 and out_score1 are scores of current module
	my $in_score1 = 0;
	my $out_score1 = 0;
	$Si = $Si/2;
	$e1 = $e1/2;
	if($e1 > 0){
		$in_score1 = $Si/$e1;
	}
	if($e2 > 0){
		$out_score1 = $So/$e2;
	}

#	my $Sd1 = $in_score1 - $out_score1; # score of current module
	my $Sd = $out_score1 - $in_score1; # score of current module

	return($Sd, $Si, $So, $e1, $e2);
}

sub cal_add_mod_score{
	my ($mod_ref, $new_node, $Si, $So, $e1, $e2) = @_;
	
	my $Sii = 0;
	my $Soo = 0;
	my $ee1 = 0;
	my $ee2 = 0;

	foreach my $g1(keys %{$FLN{$new_node}}){
		if(exists ${$mod_ref}{$g1}){
			$Sii = $Sii + $FLN{$g1}{$new_node};
			$ee1 = $ee1 + 1;
			#$So = $So - $FLN{$g1}{$new_node};
			#$e2 = $e2 - 1;
		}else{
			$Soo = $Soo + $FLN{$g1}{$new_node};
			$ee2 = $ee2 + 1;
			#$So = $So + $FLN{$g1}{$new_node};
			#$e2 = $e2 + 1;
		}
	}

	$Si = $Si+$Sii;
	$So = $So-$Sii+$Soo;
	$e1 = $e1+$ee1;
	$e2 = $e2+$ee2-$ee1;

	my $in_score1 = 0;
	my $out_score1 = 0;
	if($e1 > 0){
		$in_score1 = $Si/$e1;
		#$in_score1 = ($Si+$Sii)/($e1+$ee1);
	}
	if($e2 > 0){
		$out_score1 = $So/$e2;
		#$out_score1 = ($So-$Sii+$Soo)/($e2+$ee2-$ee1);
	}

#print "add:$in_score1\t$out_score1\n";

#	my $Sd1 = $in_score1 - $out_score1; # score of current module
	my $Sd = $out_score1 - $in_score1; # score of current module

	return($Sd, $Si, $So, $e1, $e2);
}

sub cal_toss_mod_score{
	my ($mod_ref, $toss_node, $Si, $So, $e1, $e2) = @_;
	
	my $Sii = 0;
	my $Soo = 0;
	my $ee1 = 0;
	my $ee2 = 0;
	
	foreach my $g1(keys %{$FLN{$toss_node}}){
		if(exists ${$mod_ref}{$g1}){
			$Sii = $Sii + $FLN{$g1}{$toss_node};
			$ee1 = $ee1 + 1;
=head
			$Si = $Si - $FLN{$g1}{$toss_node};
			$e1 = $e1 - 1;
			$So = $So + $FLN{$g1}{$toss_node};
			$e2 = $e2 + 1;
=cut
		}else{
			$Soo = $Soo + $FLN{$g1}{$toss_node};
			$ee2 = $ee2 + 1;
#			$So = $So - $FLN{$g1}{$toss_node};
#			$e2 = $e2 - 1;
		}
	}

	$Si = $Si - $Sii;
	$So = $So+$Sii-$Soo;
	$e1 = $e1-$ee1;
	$e2 = $e2+$ee1-$ee2;

	my $in_score1 = 0;
	my $out_score1 = 0;
	if($e1 > 0){
		#$in_score1 = ($Si-$Sii)/($e1-$ee1);
		$in_score1 = $Si/$e1;
	}
	if($e2 > 0){
		#$out_score1 = ($So+$Sii-$Soo)/($e2+$ee1-$ee2);
		$out_score1 = $So/$e2;
	}

#	my $Sd1 = $in_score1 - $out_score1; # score of current module
	my $Sd = $out_score1 - $in_score1; # score of current module

	return($Sd, $Si, $So, $e1, $e2);
}
