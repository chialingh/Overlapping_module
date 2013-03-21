#!/usr/bin/perl

#7625
for($i=1860; $i<=7625; $i++){
	system("perl merge_modules_test_number_union.pl 0.1 $i >> Module_number_test_union.txt");
}
