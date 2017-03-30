#!/usr/bin/perl
use warnings; use strict;
use Getopt::Long;
################################Get Sequences#############################################
############################################################################################
##########################################################################################
my $cut_flag=0;
my $sequences_fname;
my $output_prefix;
my $remove_Ns_flag=0;
my $z_syn_genome_cpynumber=10;
GetOptions(	"c" => \$cut_flag,
		"n" => \$remove_Ns_flag,
		"s=s" => \$sequences_fname,
		"o=s" => \$output_prefix,
		"z=i" => \$z_syn_genome_cpynumber);
# -c 1 turns on +- 1kb substring argument   -f is the cps file location -s is the sequence file location
if(!(-e $sequences_fname)){die "No Sequences file location found. Please give file location as -s argument\n";}

open(SEQ,'<',"$sequences_fname")||die "Cannot open sequences file\n";
my @file=<SEQ>;
close SEQ;
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
chomp @headers;
my $removed = 0;
#open(TRIMMED ,">","trimmed.fa")||die "Cannot open outfile \n";
for(my $i=0;$i<@sequences;++$i){
	if($cut_flag){
		############IMPORTANT, only use this loop when using 1kb+-orf
		$sequences[$i]=substr($sequences[$i],1000,length($sequences[$i])-2000);
	}
	my $problem = 0;
	# check that the protein starts with methionine and ends with a stop codon
	if(substr($sequences[$i],0,3) ne "ATG") {
		$problem = 1;	
		print "bad start " . $i . " because of start \n";	
	}
	my $end = substr($sequences[$i],(length($sequences[$i])-3),3);
	if(!(($end eq "TAA") || ($end eq "TAG") || ($end eq "TGA"))) {
		$problem = 1;	
		print "bad end " . $i . " because of " . $end . "\n";	
	}
	if((($remove_Ns_flag) && ($sequences[$i]=~m/N/)) || $problem) {
		#print "removing header " . $headers[$i] . "\n";
		$removed++;
		splice(@sequences,$i,1);
		splice(@headers,$i,1);
		--$i;
	} else {
		#print TRIMMED $headers[$i] . "\n" . $sequences[$i] . "\n";
	}
}

# read in the codon pairs and set up the maps
my $k=6;
my %cp_expected_frm;
my %cpair_counts_frm;

open(CP,"<","codon_pair.txt")||die "Cannot open  codon pair file\n";
while(my $line=<CP>){
	chomp $line;
	$cp_expected_frm{$line}=[];
	$cpair_counts_frm{$line}=[];
	for(my $i=0;$i<3;++$i){
		${$cp_expected_frm{$line}}[$i]=1;##Pseudocount of 1
		${$cpair_counts_frm{$line}}[$i]=1;
	}
}
close CP;

################################Get Actual Counts CPs###################
my $seq_with_problem = 0;
my $total_seq_with_problem = 0;
for(my $i=0;$i<@sequences;++$i){
	my $has_problem = 0;
	for(my $j=0;$j<length($sequences[$i])-$k;$j=$j+3){
		my $cod = substr($sequences[$i],$j,3);
		#if(($cod eq "TAA") || ($cod eq "TAG") || ($cod eq "TGA")) {
			#print TRIMMED "Stop Codon in the reading frame!!!!\n";
			#print TRIMMED $headers[$i] . "\t" . $j . "\t" . $cod . "\n";
			#print TRIMMED $sequences[$i] . "\n";
			#$has_problem++;
		#}
		#if(($cod eq "ATG") && ($j > 0)) {
		#	print TRIMMED "Start Codon further in reading frame!!!!\n";
		#	print TRIMMED $headers[$i] . "\t" . $j . "\t" . $cod . "\n";
		#	print TRIMMED $sequences[$i] . "\n";
		#	$has_problem++;
		#}
		if(defined($cpair_counts_frm{substr($sequences[$i],$j,$k)})) {
			++${$cpair_counts_frm{substr($sequences[$i],$j,$k)}}[0];
		}
	}
	for(my $j=1;$j<length($sequences[$i])-$k;$j=$j+3){
		if(defined($cpair_counts_frm{substr($sequences[$i],$j,$k)})) {
			++${$cpair_counts_frm{substr($sequences[$i],$j,$k)}}[1];
		}
	}
	for(my $j=2;$j<length($sequences[$i])-$k;$j=$j+3){
		if(defined($cpair_counts_frm{substr($sequences[$i],$j,$k)})) {
			++${$cpair_counts_frm{substr($sequences[$i],$j,$k)}}[2];
		}
	}
	#if($has_problem > 0) {
	#	$seq_with_problem++;
	#	$total_seq_with_problem += $has_problem;
	#	print TRIMMED $headers[$i] . "\n";
	#}
}
#print TRIMMED $seq_with_problem . " orfs have problems with stop in the middle\n"; 
#print TRIMMED $total_seq_with_problem . " total problems with stop in the middle\n"; 

##########################################Generate Randomized Genome#################
#####################################################################################
#####################################################################################
my %syn_cds;
open (CAA,'<',"./syn_codons.txt")||die "Cannot open synonomyous codons file\n";
while (my $line=<CAA>){
	chomp $line;
	my @a=split(/\t/,$line);
	my $aa=pop(@a);
	foreach(@a){
		$syn_cds{$_}=\@a;
	}
}
close CAA;
	
for(my $i=0;$i<(@sequences);++$i){#create synthetic genome which preserves amino acid structure of each sequence but randomly chooses syn codons at each position to eliminate codon pair bias
	my $l=length($sequences[$i])-3;
	my $seq=$sequences[$i];
	for(my $r=0;$r<$z_syn_genome_cpynumber;++$r){ #Generates many copies of each gene shuffled independently
		my %unique;
		my @iterate_arr=grep { ! $unique{$_} ++ } values %syn_cds;
		foreach my $syns( @iterate_arr ){
			my @index;my @permute;
			my $cmd = join '|', map { quotemeta $_ } @{$syns};
			while($sequences[$i] =~ m/\G([ACGT]{3})+?($cmd)/g){
				pos($sequences[$i])=$-[2];##resets position of next search start to appropriate location
				push(@index,$-[2]);
			}
			if (scalar(@index) <= 1){next;}
			@permute=@index;
			&fisher_yates_shuffle(\@permute);
			for(my $j=0;$j<@index;++$j){
				unless($permute[$j] eq $index[$j]){
					substr($seq,$permute[$j],3,substr($sequences[$i],$index[$j],3));
					substr($seq,$index[$j],3,substr($sequences[$i],$permute[$j],3));
				}
			}
		}	
		################################Add to Expected Counts CPs###################
		for(my $j=0;$j<length($seq)-$k;$j=$j+3){
			if(defined($cp_expected_frm{substr($seq,$j,$k)})) {
				++${$cp_expected_frm{substr($seq,$j,$k)}}[0];
			}
		}
		for(my $j=1;$j<length($seq)-$k;$j=$j+3){
			if(defined($cp_expected_frm{substr($seq,$j,$k)})) {
				++${$cp_expected_frm{substr($seq,$j,$k)}}[1];
			}
		}
		for(my $j=2;$j<length($seq)-$k;$j=$j+3){
			if(defined($cp_expected_frm{substr($seq,$j,$k)})) {
				++${$cp_expected_frm{substr($seq,$j,$k)}}[2];
			}
		}
	}
}
##divide expected counts by the copy number of the synthetic genome so that the expected counts are normalized 
foreach my $key(keys %cpair_counts_frm){
#	print "before division:\t",$key,"\t",${$cp_expected_frm{$key}}[0],"\t",${$cp_expected_frm{$key}}[1],"\t",${$cp_expected_frm{$key}}[2],"\n";
	if (${$cp_expected_frm{$key}}[0] != 1) {
		${$cp_expected_frm{$key}}[0]/=$z_syn_genome_cpynumber;
	}
	if (${$cp_expected_frm{$key}}[1] != 1) {
		${$cp_expected_frm{$key}}[1]/=$z_syn_genome_cpynumber;
	}
	if (${$cp_expected_frm{$key}}[2] != 1) {
		${$cp_expected_frm{$key}}[2]/=$z_syn_genome_cpynumber;
	}
}
my %cp_frmspc_score;
foreach my $key(keys %cpair_counts_frm){
	$cp_frmspc_score{$key}=[];##make hash of arrays with codon pair keys and three scalar entries in each array
	${$cp_frmspc_score{$key}}[0]=log(${$cpair_counts_frm{$key}}[0]/${$cp_expected_frm{$key}}[0]);
	${$cp_frmspc_score{$key}}[1]=log(${$cpair_counts_frm{$key}}[1]/${$cp_expected_frm{$key}}[1]);
	${$cp_frmspc_score{$key}}[2]=log(${$cpair_counts_frm{$key}}[2]/${$cp_expected_frm{$key}}[2]);
}


my $s="$output_prefix"."_pairwise_syn_shuffle_frm123.csv";
open( DAT ,">","$s")||die "Cannot open outfile $s\n";
foreach my $key( sort { $a cmp $b } keys %cpair_counts_frm){
	print DAT $key,",",${$cp_frmspc_score{$key}}[0],",",${$cp_frmspc_score{$key}}[1],",",${$cp_frmspc_score{$key}}[2],"\n";
}
close DAT;
my $w="$output_prefix"."_obs_exp_frm123.cts";
open (CTS,">$w")||die "Cannot open $w file\n";
foreach my $key (sort {$a cmp $b} keys %cp_expected_frm){
	print CTS $key,"\t",${$cp_expected_frm{$key}}[0],"\t",${$cp_expected_frm{$key}}[1],"\t",${$cp_expected_frm{$key}}[2],"\t",${$cpair_counts_frm{$key}}[0],"\t",${$cpair_counts_frm{$key}}[1],"\t",${$cpair_counts_frm{$key}}[2],"\n";
}
close CTS;

exit;

sub fisher_yates_shuffle {#pass array reference and the array will be randomly permuted in a random way
    my $array=shift;
        my $i; 
    for ($i=@$array;--$i;) {
                my $j=int rand ($i+1);
            next if $i==$j;
                @$array[$i,$j]=@$array[$j,$i];
        }   
}
