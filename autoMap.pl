#!/usr/bin/env perl
# ./autoMap.pl /home/someone/Zackular test.map 2
# ./autoMap.pl /home/someone/Flores,/home/someone/Zackular test3.map 2
use warnings;
use strict;
sub prefix_find;


if (@ARGV<3){
	print "automap.pl is now part of lotus2.pl, just call \"./lotus2.pl -create_map\"\n\n";
	print "Deprecated: Usage:\n./automap.pl [input dirs] [output mapping file] [paired 1/2]\n -input dirs can be one dir or several, seperated by \",\"\n";
	print " -output mapping file is relative path and name of mapping file\n -paired is the expected paired (2) or single (1) input in dirs\n\n";
	print "Example: ./automap.pl /home/singlRds1,/home/pairedRds/,home/singlRds2/ myMap.txt 1,2,1 \n\n";
	exit(0);
}
my $refDirs = $ARGV[0];
my $ofile = $ARGV[1];
my $pairedP = $ARGV[2];
my $prefix = "SMPL";
my $smpCnt = 0;
my @paired = split /,/,$pairedP;

if (@paired > 1 ){
	print "Looking for mix of paired and single input files\n";
}elsif ($paired[0]==2){
	print "Looking for paired input files\n";
} elsif ($paired[0]==1){
	print "Looking for single read input files\n";
} else {
	die "Arg 2 has to be supplied and be either \"1\" or \"2\" (single / paired read input in dirs)\n";
}
#### regex search string
my $rawFileSrchStr1 = '(R1_001|_1|\.R1)\.f[^\.]*[aq][\.gz]*$';
my $rawFileSrchStr2 = '(R2_001|_2|\.R2)\.f[^\.]*[aq][\.gz]*$';


my $rawSingSrchStr = '.f[^\.]*[aq][\.gz]*$';

my $map = "#SampleID	fastqFile\n";
my @preRDs = split(/,/,$refDirs);
if (@preRDs > 1 && @paired == 1){ #extend to full length
	for (my $i=1;$i<@preRDs;$i++){
		$paired[$i] = $paired[0];
	}
}
my ($arRDs,$pathPre) = prefix_find(\@preRDs);
my @RDs = @{$arRDs};
my $cnt=0;
foreach my $RD ( @RDs ){
	#print "DIR: $RD\n";
	opendir(DIR, "$pathPre/$RD") or die $!;
	my @pa2 =(); 
	my @pa1 = ();
	if ($paired[$cnt] == 2){
		@pa2 = sort ( grep { /$rawFileSrchStr2/ && -e "$pathPre/$RD/$_" } readdir(DIR) );	rewinddir DIR;
		@pa1 = sort ( grep { /$rawFileSrchStr1/  && -e "$pathPre/$RD/$_"} readdir(DIR) );	close(DIR);
	} else {
		@pa1 = sort ( grep { /$rawSingSrchStr/  && -e "$pathPre/$RD/$_"} readdir(DIR) );	close(DIR);
	}
	print "Found ".@pa1." samples in dir $RD\n";
	die "unequal read pairs!: ".@pa1 .",". @pa2."\n" if (@pa1 != @pa2 && $paired[$cnt] == 2);
	print "Assuming unpaired reads, as can not find signature for pair 2 files\n" if (@pa2 == 0 && @pa1 >0 && $paired[$cnt] == 2);
	#add to map
	foreach (my $i=0; $i< @pa1; $i++){
		if (@pa2 > 0){
			$map .= $prefix.$smpCnt."\t".$RD."".$pa1[$i].",".$RD."".$pa2[$i]."\n";
		} else {
			$map .= $prefix.$smpCnt."\t".$RD."".$pa1[$i]."\n";
		}
		$smpCnt ++;
	}
	$cnt++;
	
}
open O,">$ofile" or die "Can't open $ofile\n";
print O $map;
close O;

print "Map is $ofile\n";
print "Please check that all files required are present.\n";
print "==========================\nRun: ./lotus.pl -m $ofile -i $pathPre/ -s ..\n==========================\n";
exit(0);



sub prefix_find($){
	my ($ar) = @_;
	my @rds = @{$ar}; my @newRds = @rds;
	my $first = $rds[0]; 
	for (my $i=0;$i<@rds; $i++){
		my $second = $rds[$i];
		"$first\0$second" =~ m/^(.*\/).*\0\1/s;
		$first = $1;
	}
	
	print "$first common prefix\n";
	for (my $i=0; $i<@newRds; $i++){
		$newRds[$i] =~ s/^$first//;
		$newRds[$i] .= "/" if ($newRds[$i] !~ m/\/$/ && length($newRds[$i] ) != 0 );
	}
	return (\@newRds,$first);
}














