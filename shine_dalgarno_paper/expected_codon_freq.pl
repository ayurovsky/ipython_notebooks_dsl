#!/usr/bin/perl
use warnings; use strict;
##program finds expected codon pair frequency for inputted sequences
##for specific frames
my @file=<>;
my @headers;
##Findnumberofsequences
my $count=0;
foreach my $line(@file){
	$count+=()=$line=~/>/;
}
my @sequences=(('')x$count);
my $current_line;
for (my $i=0; $i<scalar @file; ++$i){#findfirstline
	if($file[$i]=~/>/){
		$current_line=$i+1;
		$i=scalar@file;
	}   
}
for(my $i=0;$i<$count;++$i){
	$headers[$i]=$file[$current_line-1];
	while($current_line<scalar @file && $file[$current_line]!~/>/){
		chomp $file[$current_line];
		$sequences[$i]="$sequences[$i]"."$file[$current_line]";
		++$current_line;
	}   
	++$current_line;
}


my $k=3;
my %ktc;
############initialize %ktc
open(CP,"<","./codons.txt")||die "Cannot open file\n";
my @cps=<CP>;
close CP;
chomp @cps;
foreach(@cps){
	$ktc{$_}=0;
	#print $_ . "\t";
}
#print "\n";
##########
my %syn_cds;
open (CAA,'<',"./syn_codons.txt")||die "Cannot open synonomyous codons file\n";
while (my $line=<CAA>){
    chomp $line;
    my @a=split(/\t/,$line);
    my $aa=pop(@a);
    foreach(@a){
        $syn_cds{$_}=\@a;
		#print $_ . "\t";
    }
}
#print "\n";
close CAA;
##########
my $ct=0;
my $frm=0;
for(my $i=0;$i<@sequences;++$i){
	for(my $j=$frm;$j<length($sequences[$i])-$frm;$j=$j+3){
		if(!defined($ktc{substr($sequences[$i],$j,$k)})) {
			# not one of the codons
			next;
		}
		++$ktc{substr($sequences[$i],$j,$k)};
		++$ct;
	}
}
#while ( ($key, $value) = each %ktc ){
#	  print "$key\t$value\n";
#}

#foreach my $key ( sort keys %ktc ){
	#  print  "$key\t$ktc{$key}\n";
#	  print  "$key\t",$ktc{$key}/$ct,"\n";
#}

my %cd_f;
foreach my $key(keys %ktc){
	$cd_f{$key}=$ktc{$key}/$ct;
#	print $key . "\t";
}
#print "\n";

foreach my $key (sort { $cd_f{$b} <=> $cd_f{$a} } keys %cd_f) {
	print $key . "\t";
	printf "%.5f\t\t(",$cd_f{$key};
	foreach(sort { $cd_f{$b} <=> $cd_f{$a} } @{$syn_cds{$key}}) {
		printf "%.5f ",$cd_f{$_};
	}
	print ")\n";
}
#print $ct; 
