#!/usr/bin/perl
use strict;

my $upstream_length = 40;
if (scalar(@ARGV) != 3) { 
    die("This script requires 3 arguments:\n\tgenes.fa\n\tgenome_reference_file.fa\n\toutput_file\n");
}
open(GENE, $ARGV[0]) || die("Cannot open the genes of interest file\n");

my $file = $ARGV[1];
if ($file =~ /.gz$/) {
        open(REF, "gunzip -c $file |") || die "can’t open genome reference $file";
}
else {
        open(REF, $file) || die "can’t open genome reference $file";
}


open(OUT, ">" . $ARGV[2]) || die("Cannot open the output file\n");


# first read in the genome reference file into a single string # for now assume a single chromosome
my $ref_chr = "nochr";
my $reference_string = "nostring";
while(my $line = <REF>) {
	chomp $line;
	if ($line =~ '^>(.*)') {
		if ($ref_chr ne "nochr") {
			#print OUT "error: more than one chromosome, fix code\n";
			#$REF_HASH{$ref_chr} = $reference_string;
			last;
		}
		$ref_chr = $1;
		$reference_string = "";
	} else {
		$reference_string = $reference_string . $line;
	}
}
# do the last line
#$REF_HASH{"ref_chr"} = $reference_string;
#print OUT $reference_string;

# print out a random name for now 
print OUT ">genes\n";

# now read the genes files, and get the substrings
while(my $line = <GENE>) {
	chomp $line;
	#if ($line =~ '^>.*:([c\d]*)-(.*) .*') { # we got the definition line
	if ($line =~ '^>.*') { # we got the definition line
		my @words = split(' ', $line);
		my $idx = substr($words[3], 10, -1);
		my $reverse = 0;
		#if ($idx =~ 'c(.*)') {
		if ($idx =~ 'complement(.*)') {
			$reverse = 1;
			$idx = $1;
			if ($idx =~ '(\d+)\.\.(\d+)') {
				$idx = $2;
			}
		}
		elsif ($idx =~ '(\d+)\.\.(\d+)') {
			$idx = $1;
		}
		# indexing is -1
		$idx -= 1;
		if ($reverse) {
		#	if (substr($reference_string,$idx-2,3) ne "CAT") {	
		#		print OUT "Problem with: " . substr($reference_string,$idx-2,3) . "\t" . $idx . "\n";
		#	} else {
			my $before_seq = substr($reference_string,$idx+1,40);
        		my $reverse_complement = reverse($before_seq); # reverse the string
        		$reverse_complement =~ tr/ACGTacgt/TGCAtgca/; # get the complement
				print OUT $reverse_complement . "|\n";
		#	}
		} else { # not reversed
		#	if (substr($reference_string,$idx,3) ne "ATG") {	
		#		print OUT "Problem with: " . substr($reference_string,$idx,3) . "\t" . $idx . "\n";
		#	} else {
				my $before_seq = substr($reference_string,$idx-40,40);
				print OUT $before_seq . "|\n";
		#	}

		}
	}
}
