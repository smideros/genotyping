#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

# USE: from a UNIX terminal:
# perl hmp2rqtl.pl sourcefile
# input tassel hapmap file
# output r/qtl read.cross csvr file
# unknown parents are changed to -, but keeps data if only one unknow parent.
# parents hard coded on columns 206 and 205. Ny001 is A

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
splice @headings, 0, 11;
unshift (@headings,'id','');
my $id = join (",",@headings);

## Create a newrow for each marker and compare to parental genotypes
my @output;
foreach (@original) {
  my @row = split("\t", $_);
  my @newrow;
  push(@newrow, $row[0],$row[2]);
  splice @row, 0, 11;
  #print scalar(@row)."\t"."$row[205]"."\t"."$row[204]\n";
  foreach my $element (@row){
    if (($row[205] eq 'N')&&($row[204]eq'N')){
      push(@newrow,'-');
    }elsif($row[205]eq'N'){
      if($element eq $row[204]){
	push(@newrow,'B');
      } elsif($element ne 'N'){
	push(@newrow,'A');
      } else {
	push(@newrow,'-');
      }
    }elsif($row[204]eq'N'){
      if($element eq $row[205]){
	push(@newrow,'A');
      } elsif($element ne 'N'){
	push(@newrow,'B');
      } else {
	push(@newrow,'-');
      }
    }
    elsif ($element eq $row[205]){
      push(@newrow, 'A');
    }
    elsif ($element eq $row[204]){
      push(@newrow, 'B');
    } else {
      push(@newrow,'-');
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
#print "$id\n";
#foreach (@output){
#  print "$_\n";
#}
## create the output file and print each of the rows from the output array
open (FILE, ">$ARGV[0].csv");
print FILE "$id\n";
foreach (@output) {
	print FILE "$_\n";
}
close FILE;
