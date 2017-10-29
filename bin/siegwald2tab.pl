#!/usr/bin/perl

use warnings;
use strict;

my $inFile = "ranking-of-metagenomics-tools-Siegwald-cleaned.tsv";

open(IN, "< $inFile");
my $lab='';
while(my $in=<IN>){
    chomp($in);
    $in =~ s/One Codex/OneCodex/;
    $in =~ s/kraken/Kraken/;

    if($in =~ /(\S+\s+level;.*error.*reads)/){
	$lab=$1;
    }
    
    
    my %method = (
	BMP     => 1,
	kraken  => 1,
	mothur  => 1,
	One     => 1,
	QIIME   => 1,
	CLARK   => 1
	);
    
    my @a=split(/\s+/, $in);
    #print "a:[[[@a]]]\n";
    if(defined($a[0]) && defined($method{$a[0]}) ){
	my @aa = split(/\t/, $in);
	#print "aa:[[[@aa]]]\n";
	my $lab2 = '';
	if($aa[0]=~/QIIME\s+(.*+)/){
	    $lab2=$1;
	    $aa[0]='QIIME';
	}
	print "Siegwald.2017\t$aa[0]\t$aa[4]\t$aa[3]\t$aa[2]\t$lab; Original values - 200(V3); Reference Database = $aa[1]; $lab2\n";
	print "Siegwald.2017\t$aa[0]\t$aa[7]\t$aa[6]\t$aa[5]\t$lab; Original values - 400(V4-V5); Reference Database = $aa[1]; $lab2\n";
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

