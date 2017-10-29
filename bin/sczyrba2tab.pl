#!/usr/bin/perl

use warnings;
use strict;

#my $inFile = shift;
my $inFile = "cami/Sczyrba-supp-tables-099127-table24.csv";
open(IN, "< $inFile");
my $lab='';
while(my $in=<IN>){
    chomp($in);
    $in =~ s/\+//;
    my %method = (
	Kraken        => 1,
	MEGAN         => 1,
	PhyloPythiaS  => 1,
	'taxator-tk'  => 1
	);
    
    my @a=split(/\s+/, $in);
    if(defined($a[0]) && defined($method{$a[0]}) ){
	my @aa = split(/\t/, $in);
	my $lab2 = $a[1];
	
	my ($sens,$ppv) = ($aa[1]/100, $aa[3]/100);
	printf "Sczyrba.2017\t$a[0]\t%0.4f\t%0.4f\t%0.4f\tnotes:$aa[7];$aa[8]; version:$lab2; in=$inFile;\n", 2*($sens*$ppv)/($sens + $ppv), $sens, $ppv;

    }
}

close(IN); 




