#!/usr/bin/perl

# Merge overlapping modules. For modules that have at least C% overlapping. 
# Using intersection nodes as seeds and searching in the union nodes for finding the representative module of two merged modules

use strict;

######
my @time = localtime(time);
print "Reading FLN:$time[2]:$time[1]:$time[0]\n";
######

my $cutoff = $ARGV[0]; # decide at what overlapping level modules should be merged. 0~1 => 0%~100%
my $theda = 0.2; # threshold when searching for the representative module

# read FLN
my $file3 = "FLN.data.human.txt"; # FLN
my %FLN;
open(IN, "/home/clhuang/lab/FLN/$file3");
while(<IN>){
	chomp $_;
	my @line = split(/\t/, $_);
	$FLN{$line[0]}{$line[1]} = abs($line[2]);
}
close IN;

######
my @time = localtime(time);
print "Reading modules:$time[2]:$time[1]:$time[0]\n";
######
my %module_list; # list of modules

my %module_gene; 
# hash that store the genes in each module. 
# The first key is the module number and the second key is the gene (Entrez ID) in each module. 

for(my $i = 1; $i <=10000; $i++){
	my $module = "/home/clhuang/lab/causal_modules/results/2012NOV06/modules/0p2/$i.txt";
	open(INF, "$module");
	while(my $g=<INF>){
		chomp $g;
		$module_gene{$i}{$g} = ();
	}
	close INF;
	$module_list{$i} = 1;
}
close IN;

######
my @time = localtime(time);
print "Strat merging modules:$time[2]:$time[1]:$time[0]\n";
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
	my @biggest;

	my @module_list_array = keys %module_list; # list of module names
	my $module_number = $#module_list_array; # number of modules - 1

######
my @time = localtime(time);
print "The $jj round.\n";
print "start:$time[2]:$time[1]:$time[0]\n";
print "Number of modules: ";
print scalar @module_list_array;
print "\n";
######


######
my @time = localtime(time);
print "Calculating overlapping:$time[2]:$time[1]:$time[0]\n";
######
	# for the first round, for each module pairs, module 1 and module 2
	# find number of common genes and percentage of common genes in each module.
	if($jj == 1){
		for(my $i = 0; $i<= $module_number - 1; $i++){
			my $m1 = $module_list_array[$i];
			for(my $j = $i + 1; $j <= $module_number; $j++){
				my $m2 = $module_list_array[$j];
				my ($value3, $value4, $value5) = find_common_genes($m1, $m2, \%module_gene);
				if($value3){
					my $pair = $m1."_".$m2; # module pair name
					$hash3{"$pair"} = $value3;
					$hash4{"$pair"} = $value4;
					$hash5{"$pair"} = $value5;
				}
			}
		}
	}else{ # for the following rounds, only calculate overlapping between new merged modules and rest of unmerged modules
		foreach my $m1(keys %merged_module_list){
			for(my $i = 0; $i<= $module_number; $i++){
				my $m2 = $module_list_array[$i];
				my ($value3, $value4, $value5) = find_common_genes($m1, $m2, \%module_gene);
				if($value3){
					my $pair = $m1."_".$m2; # module pair name
					$hash3{"$pair"} = $value3;
					$hash4{"$pair"} = $value4;
					$hash5{"$pair"} = $value5;
				}
			}
		}		
	}

######
my @time = localtime(time);
print "Merging modules:$time[2]:$time[1]:$time[0]\n";
######

	# sort module pairs according to overlapping percentage (from large to small)
	my @biggest = sort{ $hash3{$b} <=> $hash3{$a} } keys %hash3;

######
my @time = localtime(time);
my $bbb = scalar @biggest;
print "Merging $bbb modules:$time[2]:$time[1]:$time[0]\n";
######

	# finish clustering if there is no module has more than C% overlapping with another module
	last if $hash3{$biggest[0]} < $cutoff;

	# Find representative module for each needed to be merged module pair 
my $mm1 = 1;
	foreach my $pair(@biggest){
######
my @time = localtime(time);
print "Merging module $mm1:$time[2]:$time[1]:$time[0]\n";
######
		my ($module1, $module2) = split(/_/, $pair);
		# do merging only for unmerged modules
		if( exists $module_list{$module1} && exists $module_list{$module2}){ 
			my $aa = $hash4{$pair}; # percentage of common genes in module m1
			my $bb = $hash5{$pair}; # percentage of common genes in module m2
			if($aa > 2*$bb){
			# Remove module m1 if common genes are the majority of the module m1
			# which means module m1 is part of module m2
				delete $module_list{$module1};
				delete $module_gene{$module1};
			}elsif($bb > 2*$aa){
			# Remove module m2 if common genes are the majority of the module m2
			# which means module m2 is part of module m1
				delete $module_list{$module2};
				delete $module_gene{$module2};
			}else{
			# common genes are part of module m1 and module m2
			# which measn both module m1 and module m2 only capture part of the true module
			# use intersection nodes as seed and search in the union nodes for finding the true module
				my %module1_hash = %{$module_gene{$module1}}; # genes in module m1
				my %module2_hash = %{$module_gene{$module2}}; # genes in module m2
				my %intersection; # intersection nodes of module m1 and m2
				foreach(keys %module1_hash){
					$intersection{$_} = $module1_hash{$_} if exists $module2_hash{$_};
				}
				# union nodes of module m1 and m2
				my %union;
				@union{keys %module1_hash, keys %module2_hash} = ();

				# retrieve partial FLN => only contains union genes of module m1 and m2
				my %FLN1;
				foreach my $g1(keys %union){
					foreach my $g2(keys %union){
						if($FLN{$g1}{$g2}){
							$FLN1{$g1}{$g2} = $FLN{$g1}{$g2};
							$FLN1{$g2}{$g1} = $FLN{$g1}{$g2};
						}
					}
				}

				# remove module m1 and m2
				delete $module_list{$module1};
				delete $module_list{$module2};
				delete $module_gene{$module1};
				delete $module_gene{$module2};
				# new name of merged module
				my $pair_idx = $module1."X".$module2;
				$merged_module_list{$pair_idx} = ();
				my $NN = scalar keys %union; # number of union genes
				# search for representative module. %1: seeds, %2: searching field (partial FLN), %3: size of searching field
				my $merged_module_ref = find_merged_module(\%intersection, \%FLN1, $NN); 
				my %merged_module = %$merged_module_ref; # representative module
				%{$module_gene{$pair_idx}} = %merged_module; # update %module_gene with the new representative module
				undef %FLN1;
				undef %union;
				undef %intersection;
				undef %module1_hash;
				undef %module2_hash;
			}
		}
		$mm1 = $mm1 + 1;
	}
	# update module list;
	undef %module_list;
	my %module_list;
	@module_list{keys %module_gene} = 1;
	$jj = $jj + 1;
}

######
my @time = localtime(time);
print "Printing results:$time[2]:$time[1]:$time[0]\n";
######

# print results
open(OUT, ">Merged_10000_module4_0p2_$cutoff.txt");
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
	my $m1 = shift;
	my $m2 = shift;

	my $module_gene_ref = shift;
	%module_gene = %$module_gene_ref;

	my %hash1;
	my %hash2;

	my %hash1 = %{$module_gene{$m1}}; # genes in module m1
	my %hash2 = %{$module_gene{$m2}}; # genes in module m2

	my %intersection; # intersection nodes
	my %union;	  # union nodes
	my $value3 = 0;
	my $value4 = 0;
	my $value5 = 0;

	# find intersection nodes
	foreach(keys %hash1){
		$intersection{$_} = $hash1{$_} if exists $hash2{$_};
	}
	my $isn = scalar keys %intersection; # number of intersection nodes

	# find union nodes
	@union{keys %hash1, keys %hash2} = ();
	my $un = scalar keys %union; # number of union nodes

	# only store qualified module pairs
	if($isn/$un > $cutoff){
		my $s1 = scalar keys %hash1; # number of genes in module m1, size of module m1
		my $s2 = scalar keys %hash2; # number of genes in module m2, size of module m2
		$value3 = $isn/$un; # For module pair m1_m2, overlapping percentage of module pairs = intersection/union
		$value4 = $isn/$s1;  # For module pair m1_m2, percentage of common genes in module m1 = intersection/size of module m1
		$value5 = $isn/$s2;  # For module pair m1_m2, percentage of common genes in module m2 = intersection/size of module m2
	}
	return($value3, $value4, $value5);
}

sub find_merged_module{
	my $ins_ref = shift;  # reference of hash of intersection nodes
	my $fln1_ref = shift; # reference of hash of union FLN
	my $NN = shift;	      # size of union FLN

	# read intersection nodes as seeds
	my %nodes = %$ins_ref; # checked nodes
	my %module = %$ins_ref; # first module

	# partial FLN (Union FLN)
	my %FLN1 = %$fln1_ref;	

	my %new_ns; # neighbors of first layer nodes
	foreach my $g1(keys %module){
		foreach my $g2(keys %{$FLN1{$g1}}){
			unless($module{$g2}){ # if new neighbor g2 is not in current module
				$new_ns{$g2} = 1;
			}
		}
	}


	my $N = scalar keys %nodes; # number of checked nodes

	# stop if all nodes in the union FLN are checked
	while($N < $NN){
		my @all_ns = keys %new_ns;	# new candidate members

		# evaluate new neighbores
		# %1: new candidate members, %2: checked nodes, %3: current module, %4: union FLN
		my ($new_mem, $nodes) = eval_node(\@all_ns, \%nodes, \%module, \%FLN1);

		my %new_mem_hash = %$new_mem; # new member
		my %nodes = %$nodes;	      # checked nodes

		# stop if there is no new member
		last if scalar keys %new_mem_hash == 0;

		# add new members into current module
		foreach my $m(keys %new_mem_hash){
			$module{$m} = 1;
		}

		@all_ns = ();
		my %new_ns; # new candidate members, i.e. neighbors of nodes in current module
		foreach my $g1(keys %module){
			foreach my $g2(keys %{$FLN1{$g1}}){
				unless($module{$g2}){ # if new neighbor g2 is not in current module
					$new_ns{$g2} = 1;
				}
			}
		}
		last if scalar keys %new_ns == 0; # stop if there is no new neighbor
	}

	return(\%module);			
}


sub eval_node{
	my $ns_ref = shift;
	my $nd_ref = shift;
	my $mod_ref = shift;
	my $fln1_ref = shift;

	my @new_ns = @$ns_ref; # new candidate members
	my %nd = %$nd_ref; # checked nodes
	my %mod = %$mod_ref; # current module
	my %FLN1 = %$fln1_ref; # partial FLN (Union FLN)

	my %members; # new members

	#########################################
	#					#
	# Calculate score of current module: Sd #
	#					#
	#########################################
	my $Si1 = 0;
	my $So1 = 0;
	my $e1 = 0;
	my $e2 = 0;
	foreach my $g1(keys %mod){ # for every node in current module
		foreach my $g2(keys %{$FLN1{$g1}}){ # look for their neighbors on FLN
			if(exists $mod{$g2}){ # if a neighbor is inside module, calculate Si
				$Si1 = $Si1 + $FLN1{$g1}{$g2};
				$e1 = $e1 + 1;
			}else{ # if a neighbor is outside module, calculate So
				$So1 = $So1 + $FLN1{$g1}{$g2};
				$e2 = $e2 + 1;
			}
		}
		$nd{$g1} = 1;
	}

	# stop if there is no outside neighbors
	last if $e2 == 0;

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
		my $ee1 = 0;
		my $ee2 = 0;
		my $g1 = $new_ns[$i];
		$nd{$g1} = 1;
		foreach my $g2(keys %{$FLN1{$g1}}){ # for every neighbor of g1
			# examine if neighbor g2 is in the current module
			if($mod{$g2}){ 
				$Si2 = $Si2 + $FLN1{$g1}{$g2};
				$ee1 = $ee1 + 1;
			}else{
				$So2 = $So2 + $FLN1{$g1}{$g2};
				$ee2 = $ee2 + 1;
			}
		}

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
	}
	return (\%members, \%nd);
}

