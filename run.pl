#!/usr/bin/perl

for($i=1; $i<=10000; $i++){
	#system("perl find_driver_modules_NeXO.pl $i.sds 0.1");
	system("perl find_driver_modules.pl $i.sds 0.2");
}
