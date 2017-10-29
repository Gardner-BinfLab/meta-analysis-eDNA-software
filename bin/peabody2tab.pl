#!/usr/bin/perl

use warnings;
use strict;

my $inFile = "ranking-of-metagenomics-tools-Peabody.tsv";

open(IN, "< $inFile");
my $lab='';
while(my $in=<IN>){
    chomp($in);
    
    my %method = (
	CARMA3       => 1,
	CLARK        => 1,
	DiScRIBinATE => 1,
	Filtered     => 1,
	Kraken       => 1,
	MetaBin      => 1,
	metaCV       => 1,
	MetaPhyler   => 1,
	phymmBL      => 1,
	RITA         => 1,
	TACOA        => 1
	);
    
    my @a=split(/\s+/, $in);
    if(defined($a[0]) && defined($method{$a[0]}) ){
	my @aa = split(/\t/, $in);
	my $lab2 = '';
	if($aa[0]=~/Filtered/){
	    $lab2='Filtered Kraken';
	    $aa[0]='Kraken';
	}
	elsif($aa[0]=~/DiScRIBinATE\s+(.*+)/){
	    $lab2=$aa[0];
	    $aa[0]='DiScRIBinATE';
	}

	print "Peabody.2015\t$aa[0]\t$aa[5]\t$aa[3]\t$aa[2]\tcladeExclusionLevel:$aa[1]; $lab2\n";

    }
}

close(IN); 




