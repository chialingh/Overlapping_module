#!/usr/bin/perl

# Merge overlapping modules. For modules that have at least C% overlapping. 
# Using intersection nodes as seeds and searching in the union nodes for finding the representative module of two merged modules

use strict;
use warnings;
use Math::BigInt;


open(LOG, ">log_NeXO");
######
my @time = localtime(time);
print LOG "Reading FLN:$time[2]:$time[1]:$time[0]\n";
######

my $cutoff = $ARGV[0]; # decide at what overlapping level modules should be merged. 0~1 => 0%~100%
my $theda = 0; # threshold when searching for the representative module

# read FLN
my %FLN;
open(IN, "/home/clhuang/lab/databases/NeXO/nbt.pairs");
while(<IN>){
	chomp $_;
	my @line = split(/\t/, $_);
	$FLN{$line[0]}{$line[1]} = 1;
	$FLN{$line[1]}{$line[0]} = 1;
}
close IN;

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
				# $module_list{$pair_idx} = ();
				my $NN = scalar keys %union; # number of union genes

				# search for representative module. %1: seeds, %2: searching field (partial FLN), %3: size of searching field
				my $merged_module_ref = find_merged_module(\%intersection, \%union, $NN); 
				# my %merged_module = %$merged_module_ref; # representative module
				%{$module_gene{$pair_idx}} = %$merged_module_ref; # update %module_gene with the new representative module

				undef %union;
				undef %intersection;
#				undef %module1_hash;
#				undef %module2_hash;
			}
		}
		$mm1 = $mm1 + 1;
	}

	# update module list;
#	my %module_list;
#	@module_list{keys %module_gene} = 1;
	$jj = $jj + 1;
}

######
@time = localtime(time);
print LOG "Printing results:$time[2]:$time[1]:$time[0]\n";
######

close LOG;

# print results
open(OUT, ">Merged1_NeXO_0_$cutoff.txt");
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

sub find_merged_module{
	my ($ins_ref, $uni_ref, $NN) = @_;

	# read intersection nodes as seeds

	my %new_ns; # neighbors of first layer nodes
	foreach my $g1(keys %$uni_ref){
		foreach my $g2 (keys %{ $FLN{$g1} } ){
			if( ${$uni_ref}{$g2} ){
				unless( ${$ins_ref}{$g2} ){
					$new_ns{$g2} = 1;
				}
			}
		}
	}

	my $N = scalar keys %$ins_ref; # number of checked nodes

	my $nodes_ref = $ins_ref;

	# stop if all nodes in the union FLN are checked
	while($N < $NN){
		my @all_ns = keys %new_ns;	# new candidate members

		# evaluate new neighbores
		# %1: new candidate members, %2: checked nodes, %3: current module, %4: union FLN
		my ($new_mem, $nodes_ref) = eval_node(\@all_ns, $nodes_ref, $ins_ref, $uni_ref);

		#my %new_mem_hash = %$new_mem; # new member
		#my %nodes = %$nodes_ref;	      # checked nodes

		# stop if there is no new member
		last if scalar keys %$new_mem == 0;

		# add new members into current module
		foreach my $m(keys %$new_mem){
			${$ins_ref}{$m} = 1;
		}

		@all_ns = ();
		my %new_ns; # new candidate members, i.e. neighbors of nodes in current module
		foreach my $g1(keys %$ins_ref){
			foreach my $g2 (keys %{$FLN{$g1}} ){
				if( ${$uni_ref}{$g2} ){
					unless(${$ins_ref}{$g2}){ # if new neighbor g2 is not in current module
						$new_ns{$g2} = 1;
					}
				}
			}
		}
		last if scalar keys %new_ns == 0; # stop if there is no new neighbor
	}

	#return(\%module);			
	return($ins_ref);			
}


sub eval_node{
	my ($ns_ref, $nd_ref, $mod_ref, $uni_ref) = @_;

	my @new_ns = @$ns_ref; # new candidate members
	my %members; # new members

	#########################################
	#					#
	# Calculate score of current module: Sd #
	#					#
	#########################################
	my $mod_n = scalar keys %$mod_ref;
	my $new_ns_n = scalar @new_ns;
	my $Si1 = 0;
	my $So1 = 0;
	my $e1 = factorial_p($mod_n)/factorial_p($mod_n-2);;
	my $e2 = $new_ns_n * $mod_n;
	foreach my $g1(keys %$mod_ref){ # for every node in current module
		foreach my $g2(keys %{$FLN{$g1}} ){
		#foreach my $g2(keys %{${$fln1_ref}{$g1}}){ # look for their neighbors on FLN
			if( ${$uni_ref}{$g2} ){
				if(exists ${$mod_ref}{$g2}){ # if a neighbor is inside module, calculate Si
					$Si1 = $Si1 + $FLN{$g1}{$g2};
					# $Si1 = $Si1 + ${$fln1_ref}{$g1}{$g2};
					#$e1 = $e1 + 1;
				}else{ # if a neighbor is outside module, calculate So
					$So1 = $So1 + $FLN{$g1}{$g2};
					# $So1 = $So1 + ${$fln1_ref}{$g1}{$g2};
					#$e2 = $e2 + 1;
				}
			}
		}
		${$nd_ref}{$g1} = 1;
	}

	# stop if there is no outside neighbors
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
	# score of current module: Sd1
	my $Sd1 = $in_score1 - $out_score1;

	##############################################
	#					     #
	# Calculate score of module with new node g1 #
	#					     #
	##############################################
	# for every new node g1, examine if adding g1 will increase the total score Sd of module
	for(my $i = 0; $i <= $#new_ns; $i++){
		my $Si2 = 0;
		my $So2 = 0;
		my $ee1 = factorial_p($mod_n + 1)/factorial_p($mod_n + 1 - 2);
		my $g1 = $new_ns[$i];
		${$nd_ref}{$g1} = 1;
		foreach my $g2 (keys %{$FLN{$g1}} ){
		# foreach my $g2(keys %{${$fln1_ref}{$g1}}){ # for every neighbor of g1
			# examine if neighbor g2 is in the current module
			if(${$mod_ref}{$g2}){ 
				$Si2 = $Si2 + $FLN{$g1}{$g2};
				# $Si2 = $Si2 + ${$fln1_ref}{$g1}{$g2};
				# $ee1 = $ee1 + 1;
			}else{
				$So2 = $So2 + $FLN{$g1}{$g2};
				# $So2 = $So2 + ${$fln1_ref}{$g1}{$g2};
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
=head
		my $in_score2 = 0;
		my $out_score2 = 0;
		if($e1+$ee1 > 0){
			$in_score2 = ($Si1+$Si2)/($e1+$ee1);
		}
		if($e2+$ee2 > 0){
			$out_score2 = ($So1-$Si2+$So2)/($e2+$ee2);
		}
		
		# score of module with new node g1: Sd2
		my $Sd2 = $in_score2 - $out_score2;

		# add g1 to new member if its addition can increase the module score by theda
		if($Sd2 - $Sd1 > $theda){
			$members{$g1} = 1;
		}
=cut
	}
	return (\%members, $nd_ref);
}

sub factorial_p{
        my $n = shift;
        my $fact = Math::BigInt->new($n);
        return $fact->bfac();
}

