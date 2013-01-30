#!/usr/bin/perl

use strict;

my $cutoff = $ARGV[0];
my $theda = 0.2;

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


my %module_list;
my %module_gene;
for(my $i = 1; $i <=10; $i++){
	my $module = "modules/0p2/$i.txt";
	open(INF, "$module");
	while(my $g=<INF>){
		chomp $g;
		$module_gene{$i}{$g} = ();
	}
	close INF;
	$module_list{$i} = 1;
}
close IN;

my %hash3;

my @biggest;
$biggest[0] = 1;
$hash3{$biggest[0]} = 1;

while($hash3{$biggest[0]} >= $cutoff){
	my %hash1; # module1
	my %hash2; # module2
	my %hash3; # overlapping percentage of module1 and module2
	my %hash4; # percentage of common genes in module1
	my %hash5; # percentage of common genes in module2

	my @module_list_array = keys %module_list;
	my $module_number = $#module_list_array;

	for(my $i = 0; $i <= $module_number - 1; $i++){
		my $m1 = $module_list_array[$i];
		%hash1 = %{$module_gene{$m1}};
		for(my $j = $i + 1; $j <= $module_number; $j++){
			my $m2 = $module_list_array[$j];	
			%hash2 = %{$module_gene{$m2}};
			my %intersection;
			my %union;

			foreach(keys %hash1){
				$intersection{$_} = $hash1{$_} if exists $hash2{$_};
			}
			my $isn = scalar(keys %intersection);

			@union{keys %hash1, keys %hash2} = ();
			my $un = scalar keys %union;

			my $pair = $m1."_".$m2;
			$hash3{"$pair"} = $isn/$un;

			my $s1 = scalar keys %hash1;
			my $s2 = scalar keys %hash2;
			$hash4{"$pair"} = $isn/$s1;
			$hash5{"$pair"} = $isn/$s2;
		}
	}

	my @biggest = sort{ $hash3{$b} <=> $hash3{$a} } keys %hash3; # sorted module pairs according to overlapping percentage (from large to small)

	last if $hash3{$biggest[0]} < $cutoff; # finish clustering if there is no module has more than 50% overlapping with another module

	foreach my $pair(@biggest){
		if($hash3{$pair} >= $cutoff){
			my ($module1, $module2) = split(/_/, $pair);
			# do merging only for unmerged modules
			if( exists $module_list{$module1} && exists $module_list{$module2}){ 
				my $aa = $hash4{$pair};
				my $bb = $hash5{$pair};
				if($aa > 2*$bb){
					delete $module_list{$module1};
					delete $module_gene{$module1};
				}elsif($bb > 2*$aa){
					delete $module_list{$module2};
					delete $module_gene{$module2};
				}else{
					my %module1_hash = %{$module_gene{$module1}};
					my %module2_hash = %{$module_gene{$module2}};
					my %intersection;
					foreach(keys %module1_hash){
						$intersection{$_} = $module1_hash{$_} if exists $module2_hash{$_};
					}
					my %union;
					@union{keys %module1_hash, keys %module2_hash} = ();

					# retrieve partial FLN
					my %FLN1;
					foreach my $g1(keys %FLN){
						foreach my $g2(keys %FLN){
							if(exists $union{$g1} && exists $union{$g2}){
								$FLN1{$g1}{$g2} = 1;
								$FLN1{$g2}{$g1} = 1;
							}
						}
					}

					delete $module_list{$module1};
					delete $module_list{$module2};
					delete $module_gene{$module1};
					delete $module_gene{$module2};
					my $pair_idx = $module1."X".$module2;
					my $NN = scalar keys %union;
					my $merged_module_ref = find_merged_module(\%intersection, \%FLN1, $NN);
					my %merged_module = %$merged_module_ref;
					%{$module_gene{$pair_idx}} = %merged_module; # update %module_gene

#foreach my $g(keys %merged_module){
#	print "$g\t";
#}
#print "\n";
				}
			}
		}
	}
	@module_list{keys %module_gene} = 1;
}

open(OUT, ">Merged_module4_10_0p2_$cutoff.txt");
print OUT "Name\tElements\n";
foreach my $k1(keys %module_gene){
	print OUT "$k1";
	foreach my $k2(keys %{$module_gene{$k1}}){
		print OUT "\t$k2";
	}
	print OUT "\n";
}
close OUT;

sub find_merged_module{

	my $ins_ref = shift;
	my $fln1_ref = shift;
	my $NN = shift;

	# read intersection as seeds
	my %nodes = %$ins_ref; # looked nodes
	my %module = %$ins_ref; # first module

	# partial FLN
	my %FLN1 = %$fln1_ref;	

	my %new_ns; # neighbors of first layer nodes
	foreach my $g1(keys %module){
		foreach my $g2(keys %{$FLN1{$g1}}){
			unless($module{$g2}){ # if new neighbor g2 is not in current module
				$new_ns{$g2} = 1;
			}
		}
	}


	my @all_ns = keys %new_ns;

	my $N = scalar keys %nodes;

	while($N < $NN){ # 21659: total number of nodes in KeggFLN_all
		my ($new_mem, $nodes) = eval_node(\@all_ns, \%nodes, \%module, \%FLN1);
		my %new_mem_hash = %$new_mem;
		my %nodes = %$nodes;
		last if scalar keys %new_mem_hash == 0;
		foreach my $m(keys %new_mem_hash){
			$module{$m} = 1;
		}

		@all_ns = ();
		my %new_ns; # neighbors of nodes in current module
		foreach my $g1(keys %module){
			foreach my $g2(keys %{$FLN1{$g1}}){
				unless($module{$g2}){ # if new neighbor g2 is not in current module
					$new_ns{$g2} = 1;
				}
			}
		}
		last if scalar keys %new_ns == 0; # stop if there is no new neighbor
		@all_ns = keys %new_ns;		
	}

	return(\%module);			
}


sub eval_node{
	my $ns_ref = shift;
	my $nd_ref = shift;
	my $mod_ref = shift;
	my $fln1_ref = shift;

	my @new_ns = @$ns_ref; # putative new neighboring nodes
	my %nd = %$nd_ref; # looked nodes
	my %mod = %$mod_ref; # current module
	my %FLN1 = %$fln1_ref; # partial FLN

	my %members; # new members

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

	my $Sd1 = $in_score1 - $out_score1; # score of current module

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
		
		my $Sd2 = $in_score2 - $out_score2; # score of module including node g1

		if($Sd2 - $Sd1 > $theda){
			$members{$g1} = 1;
		}
	}
	return (\%members, \%nd);
}

