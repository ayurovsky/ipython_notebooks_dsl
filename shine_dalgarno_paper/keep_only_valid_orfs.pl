#!/usr/bin/perl

my $read_seq = 0;
my $discarded_seq = 0;
my $to_discard_seq = 0;

if (scalar(@ARGV) != 2) {  
    die("This script requires 2 arguments:\n\tinput fasta file\n\toutput file\n");
}

my $file = $ARGV[0];
if ($file =~ /.gz$/) {
	open(IN, "gunzip -c $file |") || die "can’t open pipe to $file";
}
	else {
	open(IN, $file) || die "can’t open $file";
}

open(OUT, ">" . $ARGV[1]) || die("Cannot open the output file\n");

my $header = "";
my $seq = "";
while (my $line = <IN>) {
	chomp $line;
	if ($line =~ />/) {
		if ($header ne "") {
			my $problem = 0;
			# this is not the first header line
			# check that the protein starts with methionine and ends with a stop codon
    		#TODO: uncomment: if((substr($seq,0,3) ne "ATG" ) && (substr($seq,0,3) ne "GTG" )) {
    		if(substr($seq,0,3) ne "ATG" ) {
	        	$problem = 1;   
    		}   
    		my $end = substr($seq,(length($seq)-3),3);
    		if(!(($end eq "TAA") || ($end eq "TAG") || ($end eq "TGA"))) {
        		$problem = 1;   
			}
			if ($problem == 0) {
				my $possible_problem = 0;
				for(my $j=0;$j<length($seq)-3;$j=$j+3){
					my $cod = substr($seq,$j,3);
					if(($cod eq "TAA") || ($cod eq "TAG") || ($cod eq "TGA")) {
						print "Problem, stop codon at position " . ($j+1) . "\n";	
						$possible_problem++;
					}
				}
				if ($possible_problem == 0) {
					print OUT $header . "\n" . $seq . "\n"; 
				} else {
					$to_discard_seq++;
				}	
			} else {
				$discarded_seq++;
			}
			$read_seq++;
			$seq = "";
		}
		$header = $line;
	} else {
		$seq .= $line;
	}
}
my $problem = 0;
#TODO: uncomment: if((substr($seq,0,3) ne "ATG" ) && (substr($seq,0,3) ne "GTG" )) {
if(substr($seq,0,3) ne "ATG") {
	$problem = 1;   
}   
my $end = substr($seq,(length($seq)-3),3);
if(!(($end eq "TAA") || ($end eq "TAG") || ($end eq "TGA"))) {
	$problem = 1;
}
if ($problem == 0) {
	my $possible_problem = 0;
	for(my $j=0;$j<length($seq)-3;$j=$j+3){
		my $cod = substr($seq,$j,3);
		if(($cod eq "TAA") || ($cod eq "TAG") || ($cod eq "TGA")) {
			print "Problem, stop codon at position " . ($j+1) . "\n";	
			$possible_problem++;
		}
	}
	if ($possible_problem == 0) {
		print OUT $header . "\n" . $seq . "\n"; 
	} else {
		$to_discard_seq++;
	}	
} else {
	$discarded_seq++;
}
$read_seq++;

print "read " . $read_seq . ", discarded due to bad end or start " . $discarded_seq . " and discarded due to stops in the middle: " . $to_discard_seq . "\n";
