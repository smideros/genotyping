#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

# USE: from a UNIX terminal:
# perl hmp2rqtl.pl sourcefile
# source file headers: rs chrom sample1 sample2
# Creates a consensus genotype column CG (for consensus genotype) from two genotypes.
# Keeps all that are identical AA, TT, CC, GG
# Valid genotypes in one sample that are missing in the other are kept
# eg. AN is changed to A and NG is changed to G
# output into a tab delimited unix file named <inputfile>.out:
# rs chrom sample1 sample2 CG

## Verify and process the source file
if ($#ARGV<0) {
	print "Error: no input file \n";
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
#push(@headings, "class","CG");
#my $head = join("\t", @headings);
## Create a newrow for each marker and compare to parental genotypes
my @output;
foreach (@original) {
  my @row = split("\t", $_);
  my @newrow;
  push(@newrow, $row[0],$row[2]);
  splice @row, 0, 11;
  print scalar(@row)."\t"."$row[205]"."\t"."$row[204]\n";
  foreach my $element (@row){
    #print $element."$row[205]\n";
    if ($element eq $row[205]){
      push(@newrow, 'A');
    }
    elsif ($element eq $row[204]){
      push(@newrow, 'B');
    }
    else {push(@newrow,'-');
    }
  }
  my $marker = join(",",@newrow);
  push(@output,$marker);
}
## Debuger
#print "@output\n";
## add headers to output array
#unshift(@output, $head);
## create screen output
foreach (@output){
  print "$_\n";
}
## create the output file and print each of the rows from the output array
#open (FILE, ">$ARGV[0].out");
#foreach (@output) {
#	print FILE "$_\n";
#}
#close FILE;
