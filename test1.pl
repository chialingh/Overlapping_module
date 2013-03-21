#!/usr/bin/perl

@a = (1,2,3,4,5);
@b = (1,2,3,5,6,7);

foreach $aa(@a){
	foreach $bb(@b){
		print "$aa\t$bb\n";
		if($aa == $bb){
			last;
		}else{
			print "hh:$aa\n";
		}
	}
}
