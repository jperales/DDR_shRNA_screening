#!/usr/bin/perl -W

## Author: Perales-Paton, J
## Description: 

use strict;
use FindBin qw($Bin);

## My version
my $version="0.1";

# USAGE :
# script.pl <shRNA_reference.fa> <alignment.sam>
#

# Sanity check
if(scalar(@ARGV)!=2) {
	print STDERR "Error : 2 arguments are mandatory.","\n";
	exit 1;
}

#
my $shRNA_ref=$ARGV[0];
my $SAMfl=$ARGV[1];

# Second round of sanity check
unless(-e $shRNA_ref) {
	print STDERR "Error : the shRNA_reference.fa file  does _NOT_ exist.","\n";
	exit 1;
}
unless(-e $SAMfl) {
	print STDERR "Error : The alignment SAM file does _NOT_ exist.","\n";
	exit 1;
}

# About the shRNA reference lib
open(FA,"$shRNA_ref") or die("Cannot open $shRNA_ref : $!");

my %hairpin_LIB;
while(my $line=<FA>) {
	chomp($line);
	if($line =~ /^>(.*)$/) {
		$hairpin_LIB{$1}=0;
	}
}
close(FA);




# About the SAM file
my $not_aligned=0;
open(SA,"$SAMfl") or die("Cannot open SAM file $SAMfl : $!");
while(my $line=<SA>){
	if($line =~ /^@/) {next;}
	chomp($line);
	my @fields=split("\t",$line);
	if($fields[1]==4 or $fields[2] eq '*') {
		$not_aligned++;
	} else {
		if(exists($hairpin_LIB{$fields[2]})){
			$hairpin_LIB{$fields[2]}++;
		} else {
			print STDERR "Warn : $fields[0] aligned to an unexpected contig.","\n";
			exit 1;
		}
	}
}
close(SA);
# Print the list
my @sorted_keys=sort(keys(%hairpin_LIB));
foreach my $key (@sorted_keys){
	print STDOUT $key,"\t",$hairpin_LIB{$key},"\n";

}
print STDOUT "NOTaligned","\t",$not_aligned,"\n";
exit 0;
