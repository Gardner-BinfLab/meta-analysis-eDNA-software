#!/usr/bin/perl

use warnings;
use strict;

my $inFile = 'data/bazinet/ranking-of-metagenomics-tools-Bazinet.tsv';

my %method = (
    CARMA3       => 1,
    MEGAN        => 1,
    MetaPhyler   => 1,
    'MG-RAST'    => 1,
    MLTreeMap    => 1,
    NBC          => 1,
    PhyloPythiaS => 1,
    phymmBL      => 1,
    RAIphy       => 1,
    Treephyler   => 1
    );

print "Paper\tMethod\tF1.measure\tSensitivity\tPPV\tNotes\n";
open(IN, "< $inFile");
my $lab='';
while(my $in=<IN>){
    chomp($in);

    $in =~ s/CARMA/CARMA3/g;
    $in =~ s/PhylopythiaS/PhyloPythiaS/g;
    $in =~ s/PhymmBL/phymmBL/g;
    $in =~ s/Ralphy/RAIphy/g;
    $in =~ s/RAlphy/RAIphy/g;
    
    if($in =~ /(.*dataset)/){
	$lab=$1;
    }
    
    my @a=split(/\s+/, $in);
    my @aa = split(/\t/, $in);
    #print "a:[[[@a]]]\n";
    if(defined($a[0]) && defined($method{$a[0]}) && defined($aa[3]) && defined($aa[2]) && isNumeric($aa[3]) && isNumeric($aa[2]) ){
	#print "aa:[[[@aa]]]\n";
	my($sens, $ppv) = ($aa[2]/100, $aa[3]/100);
	printf "Bazinet.2012\t$aa[0]\t%0.4f\t%0.4f\t%0.4f\t$lab; \n", (2*$sens*$ppv)/($sens+$ppv), $sens, $ppv;

    }
    
    
}




close(IN); 




######################################################################
sub isNumeric {
    my $num = shift;
    if ($num=~/^-?\d+\.?\d*$/) { 
        return 1; 
    }
    else {
        return 0;
    }
}
