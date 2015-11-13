#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

# USE: from a UNIX terminal:
# perl consg.pl sourcefile
# source file headers: rs chrom sample1 sample2
# Creates a consensus genotype column CG (for consensus genotype) from two genotypes.
# Keeps all that are identical AA, TT, CC, GG
# Valid genotypes in one sample that are missing in the other are kept
# eg. AN is changed to A and NG is changed to G
# output into a tab delimited unix file named consg.out:
# rs chrom sample1 sample2 CG
# warning: will replace another consg.out file in the folder

## Verify and process the source file
if ($#ARGV<0) {
	print "Error: no input file $! \nInput headers: rs chrm sample1 sample2 \n";
	exit;	
}
open IN1, $ARGV[0];
my @original = <IN1>;
chomp(@original);
## Debuger
#print "@original\n";
## Remove and process the header line
my $header = shift(@original);
my @headings = split("\t", $header);
push(@headings, "class","CG");
my $head = join("\t", @headings);
## Create the consensus for each SNP in the firstrow array, then merge
## the array and add it to the output array
my @output;
foreach (@original) {
  my @firstrow = split("\t", $_);
  if ($firstrow[2] eq $firstrow[3]) {
    if (($firstrow[2] eq 'T')||($firstrow[2] eq 'A')||($firstrow[2] eq 'C')||($firstrow[2] eq 'G')) {
      push(@firstrow, 'ok');
      push(@firstrow, $firstrow[2]);
    }
    else {
      push(@firstrow, 'oknATCG', $firstrow[2]);
  }}
  elsif (($firstrow[2] eq 'A') and ($firstrow[3] eq 'N')) {
    push(@firstrow, 'miss');
    push(@firstrow, $firstrow[2]);
  }
  elsif (($firstrow[2] eq 'T') and ($firstrow[3] eq 'N')) {
    push(@firstrow, 'miss');
    push(@firstrow, $firstrow[2]);
  }
  elsif (($firstrow[2] eq 'C') and ($firstrow[3] eq 'N')) {
    push(@firstrow, 'miss');
    push(@firstrow, $firstrow[2]);
  }
  elsif (($firstrow[2] eq 'G') and ($firstrow[3] eq 'N')) {
    push(@firstrow, 'miss');
    push(@firstrow, $firstrow[2]);
  }
  elsif (($firstrow[2] eq 'N') and ($firstrow[3] eq 'A'))	{
    push(@firstrow, 'miss');
    push(@firstrow, $firstrow[3]);
  }
  elsif (($firstrow[2] eq 'N') and ($firstrow[3] eq 'T'))	{
    push(@firstrow, 'miss');
    push(@firstrow, $firstrow[3]);
  }
  elsif (($firstrow[2] eq 'N') and ($firstrow[3] eq 'C'))	{
    push(@firstrow, 'miss');
    push(@firstrow, $firstrow[3]);
  }
  elsif (($firstrow[2] eq 'N') and ($firstrow[3] eq 'G'))	{
    push(@firstrow, 'miss');
    push(@firstrow, $firstrow[3]);
  }
  else {
    push(@firstrow, 'err','NA');
  }
  my $SNP = join("\t", @firstrow);
  push(@output, $SNP);
      }
## Debuger
#print "@output\n";
## add headers to output array
unshift(@output, $head);
## create screen output
foreach (@output){
  print "$_\n";
}
## create the output file and print each of the rows from the output array
open (FILE, ">consg.out");
foreach (@output) {
	print FILE "$_\n";
}
close FILE;
