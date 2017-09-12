#!/usr/bin/perl

use warnings;
use strict;

my $inFile = "data/mcintyre/mcintyre2017-supp-table.txt";

my @datasets = qw(
HC1
HC2
LC1
LC2
LC3
LC4
LC5
LC6
LC7
LC8
simHC
simMC
simLC
Buc12
CParMed48
Gut20
Hous31
Hous21
NYCSM20
Soi50
simBA525
);

my @methods = qw(CLARK CLARK-S Kraken LMAT MEGAN DMEGAN NBC Meta-ClassifierP Meta-ClassifierR);

my %inclMethods = (
    CLARK      => 1,
    'CLARK-S'  => 1,
    Kraken     => 1,
    LMAT       => 1,
    MEGAN      => 1,
    NBC        => 1
);

my (%precisionMethods, %recallMethods);

open(IN, "< $inFile");
my $lab='';
my ($methCnt, $datasetCnt)=(0,0);
while(my $in=<IN>){
    chomp($in);
    
    my @a=split(/\s+/, $in);
    
    if($in =~ /^Precision/){
	for(my $i=1; $i<scalar(@a); $i++){
	    $precisionMethods{ $methods[$i-1] } = $a[$i] if defined($methods[$i-1]);
	}
    }
    elsif($in =~ /^Recall/){
	for(my $i=1; $i<scalar(@a); $i++){
	    $recallMethods{ $methods[$i-1] } = $a[$i] if defined($methods[$i-1]);
	}
	
	foreach my $m (sort keys %inclMethods){
	    my ($sens,$ppv) = ($recallMethods{$m}/100, $precisionMethods{$m}/100);
	    printf "McIntyre.2017\t$m\t%0.4f\t%0.4f\t%0.4f\tTable3:$datasets[$datasetCnt];Precision:$precisionMethods{$m};Recall:$recallMethods{$m};\n", 2*($sens*$ppv)/($sens + $ppv), $sens, $ppv if (defined($precisionMethods{$m}) and defined($recallMethods{$m}));
	}
	$datasetCnt++;
    }
    
}

close(IN); 




