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
# output into a tab delimited unix file named <inputfile>.out:
# rs chrom sample1 sample2 CG

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
    elsif ($firstrow[2] eq 'N'){
      push(@firstrow, 'miss2', $firstrow[2]);
    }
    else {push(@firstrow, 'okhet', $firstrow[2]);
    }
  }
  elsif($firstrow[2] eq 'N') {
    if(($firstrow[3] eq 'T')||($firstrow[3] eq 'A')||($firstrow[3] eq 'C')||($firstrow[3] eq 'G')){
      push (@firstrow, 'miss', $firstrow[3]);
    }
    elsif(($firstrow[3] eq 'R')||($firstrow[3] eq 'Y')||($firstrow[3] eq 'S')||($firstrow[3] eq 'W')||($firstrow[3] eq 'M')||($firstrow[3] eq 'K')){
      push(@firstrow, 'missh', $firstrow[3]);
    }
    else{push(@firstrow, 'cerr', '?NA');
       }
  }
  elsif($firstrow[3] eq 'N') {
    if(($firstrow[2] eq 'T')||($firstrow[2] eq 'A')||($firstrow[2] eq 'C')||($firstrow[2] eq 'G')){
      push (@firstrow, 'miss', $firstrow[2]);
    }
    elsif(($firstrow[2] eq 'R')||($firstrow[2] eq 'Y')||($firstrow[2] eq 'S')||($firstrow[2] eq 'W')||($firstrow[2] eq 'M')||($firstrow[2] eq 'K')){
      push(@firstrow, 'missh', $firstrow[2]);
    }
    else{push(@firstrow, 'cerr', '?NA');
       }
  }
  elsif(($firstrow[2] eq 'R')||($firstrow[2] eq 'Y')||($firstrow[2] eq 'S')||($firstrow[2] eq 'W')||($firstrow[2] eq 'M')||($firstrow[2] eq 'K')||($firstrow[3] eq 'R')||($firstrow[3] eq 'Y')||($firstrow[3] eq 'S')||($firstrow[3] eq 'W')||($firstrow[3] eq 'M')||($firstrow[3] eq 'K')){
    push(@firstrow, 'errh', 'N');
  }
  else{push(@firstrow, 'err', 'N');
     }
  my $SNP = join("\t", @firstrow);
  push(@output, $SNP);
}
## Debuger
#print "@output\n";
## add headers to output array
unshift(@output, $head);
## create screen output
#foreach (@output){
#  print "$_\n";
#}
## create the output file and print each of the rows from the output array
open (FILE, ">$ARGV[0].out");
foreach (@output) {
	print FILE "$_\n";
}
close FILE;
## Count class ocurrences and output to screen
my @class;
foreach (@output){
  my @outrow = split("\t", $_);
  push(@class, $outrow[4]);
}
shift (@class);
my $class2 = join("\t", @class);
my $okcount = () = $class2 =~ /\bok\b/g;
my $miss2 =()= $class2 =~ /\bmiss2\b/g;
my $okhet = () = $class2 =~ /\bokhet\b/g;
my $miss = ()= $class2 =~ /\bmiss\b/g;
my $missh = () = $class2 =~ /\bmissh\b/g;
my $err = () = $class2 =~ /\berr\b/g;
my $errh = ()= $class2 =~ /\berrh\b/g;
my $cerr =()= $class2 =~ /\bcerr\b/g;
print "Counts by Class for $ARGV[0]\n";
print "ok\t" . $okcount . "\n";
print "miss2\t" . $miss2 . "\n";
print "okhet\t" . $okhet . "\n";
print "miss\t" . $miss . "\n";
print "missh\t" . $missh . "\n";
print "cerr\t" . $cerr . "\n";
print "errh\t" . $errh . "\n";
print "err\t" . $err . "\n";
