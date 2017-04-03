#!/usr/bin/perl
use warnings; use strict;
use Getopt::Long;
use List::Util qw[min max];

my $sd_seq = "";
my $cps_file = "";
GetOptions("s=s" => \$sd_seq,
			"c=s" => \$cps_file);
open(CPS,$cps_file)||die "Can not open the file with codon pair scores\n";

# get the match strings
my %matches;
if (length($sd_seq) == 6) {
	# 6/6 match
	$matches{$sd_seq} = 100;
	foreach my $pos (0..5) {
		foreach my $letter ("A","T","C","G") {
			my $seq = substr($sd_seq,0,$pos) . $letter . substr($sd_seq,$pos+1,5-$pos);	
			if (!defined($matches{$seq})) {
				$matches{$seq} = 100;
			}
		}	
	}
} elsif (length($sd_seq) == 5) {
	foreach my $letter ("A","T","C","G") {
		my $seq = $letter . $sd_seq;	
		if (!defined($matches{$seq})) {
			$matches{$seq} = 100;
		}
		$seq = $sd_seq . $letter;	
		if (!defined($matches{$seq})) {
			$matches{$seq} = 100;
		}
	}	
}
#print "codon pair,f1 cps,f2 cps,f3cps,signal score,gap score\n";
while (my $line=<CPS>) {
	chomp $line;
	my @ts = split(",",$line);	
	if (defined($matches{$ts[0]})) {
		print $ts[0] . "," . $ts[1] . "," . $ts[2] . "," . $ts[3] . "," . (($ts[1]+$ts[2]+$ts[3])/3) . "," . ($ts[1]-(($ts[2]+$ts[3])/2)) ."\n";
	}
}
