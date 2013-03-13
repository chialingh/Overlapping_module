#!/usr/bin/perl

# Merge overlapping modules. For modules that have at least C% overlapping. 
# Using intersection nodes as seeds and searching in the union nodes for finding the representative module of two merged modules

use strict;
use warnings;


open(LOG, ">log_NeXO");
######
my @time = localtime(time);
print LOG "Reading FLN:$time[2]:$time[1]:$time[0]\n";
######

my $cutoff = $ARGV[0]; # decide at what overlapping level modules should be merged. 0~1 => 0%~100%
my $theda = 0; # threshold when searching for the representative module


######
@time = localtime(time);
print LOG "Reading modules:$time[2]:$time[1]:$time[0]\n";
######
# my %module_list; # list of modules

my %module_gene; 
# hash that store the genes in each module. 
# The first key is the module number and the second key is the gene (Entrez ID) in each module. 

open(IN, "/home/clhuang/lab/causal_modules/results/Final/Overlapping_module/NeXO/modules/0/list.txt");
while(my $ff=<IN>){
	my $module = "/home/clhuang/lab/causal_modules/results/Final/Overlapping_module/NeXO/modules/0/$ff";
	$ff =~ s/\.txt//;
	chomp $ff;
	open(INF, "$module");
	while(my $g=<INF>){
		chomp $g;
		$module_gene{$ff}{$g} = ();
	}
	close INF;
}
close IN;

######
@time = localtime(time);
print LOG "Strat merging modules:$time[2]:$time[1]:$time[0]\n";
######

my %hash3;
my %merged_module_list; # merged modules

my @biggest;
$biggest[0] = 1;
$hash3{$biggest[0]} = 1;

my $jj = 1;
while($hash3{$biggest[0]} >= $cutoff){
	my %hash3; # overlapping percentage of module pairs
	my %hash4; # percentage of common genes in module1
	my %hash5; # percentage of common genes in module2
	my @biggest = ();

	my @module_list_array = keys %module_gene; # list of module names
	my $module_number = $#module_list_array; # number of modules - 1

################################################################
=head
	foreach(@module_list_array){
		print "$_\t";
	}
	print "\n";
=cut
	@time = localtime(time);
	print LOG "The $jj round.\n";
	print LOG "start:$time[2]:$time[1]:$time[0]\n";
	print LOG "Number of modules: ";
	print LOG scalar @module_list_array;
	print LOG "\n";

	@time = localtime(time);
	my $hr1 = $time[2];
	my $min1 = $time[1];
	my $sec1 = $time[0];
	print LOG "Calculating overlapping:$time[2]:$time[1]:$time[0]\n";
##################################################################

	# for the first round, for each module pairs, module 1 and module 2
	# find number of common genes and percentage of common genes in each module.
	if($jj == 1){
		for(my $i = 0; $i<= $module_number - 1; $i++){
			my $m1 = $module_list_array[$i];
			my $value3 = 0;
			my $value4 = 0;
			my $value5 = 0;
			my $pair;
			for(my $j = $i + 1; $j <= $module_number; $j++){
				my $m2 = $module_list_array[$j];
				my ($Nvalue3, $Nvalue4, $Nvalue5) = find_common_genes($m1, $m2, \%module_gene, $cutoff);
				if($Nvalue3 > $value3){
					$pair = $m1."_".$m2; # module pair name
					$value3 = $Nvalue3;
					$value4 = $Nvalue4;
					$value5 = $Nvalue5;
				}
			}
			if($pair){
				$hash3{"$pair"} = $value3;
				$hash4{"$pair"} = $value4;
				$hash5{"$pair"} = $value5;
			}
		}
	}else{ # for the following rounds, only calculate overlapping between new merged modules and rest of unmerged modules
		foreach my $m1(keys %merged_module_list){
			my $value3 = 0;
			my $value4 = 0;
			my $value5 = 0;
			my $pair;
			for(my $i = 0; $i<= $module_number; $i++){
				my $m2 = $module_list_array[$i];
				if($m1 ne $m2){
					my ($Nvalue3, $Nvalue4, $Nvalue5) = find_common_genes($m1, $m2, \%module_gene, $cutoff);
					if($Nvalue3 > $value3){
						$pair = $m1."_".$m2; # module pair name
						$value3 = $Nvalue3;
						$value4 = $Nvalue4;
						$value5 = $Nvalue5;
					}
				}
			}
#print LOG "$pair\t$value3\n";
			if($pair){
				$hash3{"$pair"} = $value3;
				$hash4{"$pair"} = $value4;
				$hash5{"$pair"} = $value5;
			}
		}		
	}

	# finish clustering if there is no module has more than C% overlapping with another module
	last if scalar keys %hash3 < 1;

	# sort module pairs according to overlapping percentage (from large to small)
	@biggest = sort{ $hash3{$b} <=> $hash3{$a} } keys %hash3;

################################################################
	@time = localtime(time);
	my $hr2 = $time[2];
	my $min2 = $time[1];
	my $sec2 = $time[0];
	my $elapse_time = 3600*($hr2-$hr1) + 60*($min2-$min1) + ($sec2-$sec1);
	print LOG "\ntime for calculating common genes: $elapse_time\n";
	my $bbb = scalar @biggest;
	print LOG "Merging $bbb modules:$time[2]:$time[1]:$time[0]\n";
	print LOG "biggest module pair: $biggest[0]\n";
	print LOG "biggest overlapping:$hash3{$biggest[0]}\n";
#################################################################

	# Find representative module for each needed to be merged module pair 

	my $mm1 = 1;

	foreach my $pair(@biggest){
		my ($module1, $module2) = split(/_/, $pair);
		# do merging only for unmerged modules
		if( exists $module_gene{$module1} && exists $module_gene{$module2}){ 
#################################################################
	@time = localtime(time);
	print LOG "Merging module $mm1-$pair: $time[2]:$time[1]:$time[0]\n";
################################################################
			my $aa = $hash4{$pair}; # percentage of common genes in module m1
			my $bb = $hash5{$pair}; # percentage of common genes in module m2
			if($aa > 2*$bb){
			# Remove module m1 if common genes are the majority of the module m1
			# which means module m1 is part of module m2
				# delete $module_list{$module1};
				delete $module_gene{$module1};
				delete $merged_module_list{$module1} if exists $merged_module_list{$module1};
			}elsif($bb > 2*$aa){
			# Remove module m2 if common genes are the majority of the module m2
			# which means module m2 is part of module m1
				# delete $module_list{$module2};
				delete $module_gene{$module2};
				delete $merged_module_list{$module2} if exists $merged_module_list{$module2};
			}else{
			# common genes are part of module m1 and module m2
			# which measn both module m1 and module m2 only capture part of the true module
			# use intersection nodes as seed and search in the union nodes for finding the true module
			#	my %module1_hash = %{$module_gene{$module1}}; # genes in module m1
			#	my %module2_hash = %{$module_gene{$module2}}; # genes in module m2

				# intersection nodes of module m1 and m2
				my %intersection; 
				foreach(keys %{$module_gene{$module1}}){
					$intersection{$_} = ${$module_gene{$module1}}{$_} if exists ${$module_gene{$module2}}{$_};
				}

				# union nodes of module m1 and m2
				my %union;
				@union{keys %{$module_gene{$module1}}, keys %{$module_gene{$module2}} } = ();

				# remove module m1 and m2
				# delete $module_list{$module1};
				# delete $module_list{$module2};
				delete $module_gene{$module1};
				delete $module_gene{$module2};
				delete $merged_module_list{$module1} if exists $merged_module_list{$module1};
				delete $merged_module_list{$module2} if exists $merged_module_list{$module2};
				# new name of merged module
				my $pair_idx = $module1."X".$module2;
print LOG "$pair_idx\n";
				$merged_module_list{$pair_idx} = ();
				%{$module_gene{$pair_idx}} = %union; # update %module_gene with the new representative module
				undef %union;
				undef %intersection;
			}
		}
		$mm1 = $mm1 + 1;
	}
	$jj = $jj + 1;
}

######
@time = localtime(time);
print LOG "Printing results:$time[2]:$time[1]:$time[0]\n";
######

close LOG;

# print results
open(OUT, ">Merged_union_0p2_$cutoff.txt");
#open(OUT, ">test123.txt");
print OUT "Name\tElements\n";
foreach my $k1(keys %module_gene){
	print OUT "$k1";
	foreach my $k2(keys %{$module_gene{$k1}}){
		print OUT "\t$k2";
	}
	print OUT "\n";
}
close OUT;

sub find_common_genes{
	my ($m1, $m2, $module_gene_ref, $cutoff) = @_;
	# %module_gene = %$module_gene_ref;

	my %intersection; # intersection nodes
	my %union;	  # union nodes
	my $value3 = 0;
	my $value4 = 0;
	my $value5 = 0;

	# find intersection nodes
	foreach(keys %{${$module_gene_ref}{$m1}}){
		$intersection{$_} = ${$module_gene_ref}{$m1}{$_} if exists ${$module_gene_ref}{$m2}{$_};
	}
	my $isn = scalar keys %intersection; # number of intersection nodes

	# find union nodes
	@union{keys %{${$module_gene_ref}{$m1}}, keys %{${$module_gene_ref}{$m2}}} = ();
	my $un = scalar keys %union; # number of union nodes

	# only store qualified module pairs
	if($isn/$un > $cutoff){
		my $s1 = scalar keys %{${$module_gene_ref}{$m1}}; # number of genes in module m1, size of module m1
		my $s2 = scalar keys %{${$module_gene_ref}{$m2}}; # number of genes in module m2, size of module m2
		$value3 = $isn/$un; # For module pair m1_m2, overlapping percentage of module pairs = intersection/union
		$value4 = $isn/$s1;  # For module pair m1_m2, percentage of common genes in module m1 = intersection/size of module m1
		$value5 = $isn/$s2;  # For module pair m1_m2, percentage of common genes in module m2 = intersection/size of module m2
	}
	return($value3, $value4, $value5);
}
