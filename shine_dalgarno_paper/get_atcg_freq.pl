#!/usr/bin/perl
use warnings; use strict;
use Getopt::Long;

my $in_file = "";
my $SD = "";
GetOptions("s=s" => \$SD,
		   "i=s" => \$in_file);
open(IN,$in_file) || die("Can't open the input file");

my $u_count = 0;
my $a_count = 0;
my $c_count = 0;
my $g_count = 0;
my $nucleotide_count = 0;
my $SD_count = 0;

while(my $line = <IN>) {
	my $search_string = <IN>;
	chomp $search_string;
	my $string_index = index($search_string,$SD,0);
	while($string_index != -1){
		$SD_count++;
		$string_index = index($search_string,$SD,$string_index+1);
	}
	$string_index = index($search_string,"T",0);
	while($string_index != -1){
		$u_count++;
		$string_index = index($search_string,"T",$string_index+1);
	}
	$string_index = index($search_string,"t",0);
	while($string_index != -1){
		$u_count++;
		$string_index = index($search_string,"t",$string_index+1);
	}
	$string_index = index($search_string,"A",0);
	while($string_index != -1){
		$a_count++;
		$string_index = index($search_string,"A",$string_index+1);
	}
	$string_index = index($search_string,"a",0);
	while($string_index != -1){
		$a_count++;
		$string_index = index($search_string,"a",$string_index+1);
	}
	$string_index = index($search_string,"C",0);
	while($string_index != -1){
		$c_count++;
		$string_index = index($search_string,"C",$string_index+1);
	}
	$string_index = index($search_string,"c",0);
	while($string_index != -1){
		$c_count++;
		$string_index = index($search_string,"c",$string_index+1);
	}
	$string_index = index($search_string,"G",0);
	while($string_index != -1){
		$g_count++;
		$string_index = index($search_string,"G",$string_index+1);
	}
	$string_index = index($search_string,"g",0);
	while($string_index != -1){
		$g_count++;
		$string_index = index($search_string,"g",$string_index+1);
	}
	$nucleotide_count += length($search_string);
}
#print "u count =" . $u_count . " and a count =" . $a_count . "\n";
#print "c count =" . $c_count . " and g count =" . $g_count . "\n";
#print "nucleotide count = " . $nucleotide_count . "\n";
my $u_frequency = $u_count/$nucleotide_count;
my $a_frequency = $a_count/$nucleotide_count;
my $c_frequency = $c_count/$nucleotide_count;
my $g_frequency = $g_count/$nucleotide_count;
#print "u frequency = " . $u_frequency . " a frequency = " . $a_frequency . "\n";
#print "c frequency = " . $c_frequency . " g frequency = " . $g_frequency . "\n";
#print "observed = " . $SD_count . "\n";

my $expected = $nucleotide_count;
for my $c (split //, $SD) {
	if ($c eq "A") {
		$expected *= $a_frequency;
	} elsif ($c eq "C") {
		$expected *= $c_frequency;
	} elsif ($c eq "G") {
		$expected *= $g_frequency;
	} else {
		$expected *= $u_frequency;
	}
}
#print "expected = " . $expected . "\n";
#print "score is = " . log($SD_count/$expected) . "\n"; 
print log($SD_count/$expected); 

