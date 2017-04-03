#!/usr/bin/perl
use warnings; use strict;
use Getopt::Long;

my $freq_file = "";
my $sd_seq = "";
GetOptions("s=s" => \$sd_seq,
			"f=s" => \$freq_file);
open(FREQ,$freq_file)||die "Can not open the file with frequencies scores\n";

my %FREQ;
while(my $line = <FREQ>){
	chomp $line;
	my @ts = split("\t",$line);
	$FREQ{$ts[0]} = $ts[1]*1000;	
}

# get real stats
my $r_product = 1;
my $r_sum = 0;
for(my $i = 0; $i <= (length($sd_seq) - 3); $i++) {
	my $codon = substr($sd_seq,$i,3);
	$r_product *= $FREQ{$codon};	
	$r_sum += $FREQ{$codon};
}
my $r_avg = $r_sum / (length($sd_seq)-2);
#print $r_product . "\t" . $r_sum . "\t" . $r_avg . "\n";

# generate all strings
my @old_list = ("A","T","C","G");
for(my $i = 1; $i < length($sd_seq); $i++) {
	my @new_list;
	foreach my $partial (@old_list) {
		push (@new_list, $partial."A");
		push (@new_list, $partial."C");
		push (@new_list, $partial."G");
		push (@new_list, $partial."T");
	}
	@old_list = @new_list;
}

# get the p_value
my $product_smaller = 0;
my $avg_smaller = 0;
my $sum_smaller = 0;
foreach my $full_string (@old_list) {
	my $product = 1;
	my $sum = 0;
	#print $full_string . "\t";
	for(my $i = 0; $i <= (length($sd_seq) - 3); $i++) {
		my $codon = substr($full_string,$i,3);
		#print $FREQ{$codon} . "\t";
		$product *= $FREQ{$codon};	
		$sum += $FREQ{$codon};
	}
	my $avg = $sum / (length($sd_seq)-2);
	#print $product . "\t" . $avg . "\n";
	if ($product <= $r_product) {
		$product_smaller++;
	}
	if ($avg <= $r_avg) {
		$avg_smaller++;
	}
	if ($sum <= $r_sum) {
		$sum_smaller++;
	}
}
my $total_num = scalar (@old_list);
#print $r_product . "\t" . $r_avg . "\n";
#print $product_smaller . "\t" . $avg_smaller . "\n";
#print "product: " . $product_smaller/$total_num . "\tavg: " . $avg_smaller/$total_num . "\n";
print $avg_smaller/$total_num;
