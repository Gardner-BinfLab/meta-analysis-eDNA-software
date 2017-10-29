#!/usr/bin/perl

use warnings;
use strict;

#my $inFile = shift;

open(IN, "< lindgreen2016-data.tsv");
my @in = <IN>;
my @meth = split(/\t/, $in[0]); 
my @sen  = split(/\t/, $in[1]); 
my @ppv  = split(/\t/, $in[2]); 
my $lab='';

    my %method = (
	Genometa    =>1,
	Kraken	    =>1,
	LMAT	    =>1,
	MetaPhlAn   =>1,
	'MG-RAST'   =>1,
	MEGAN	    =>1,
	Qiime	    =>1,
	mOTU	    =>1,
	MetaPhyler  =>1,
	'taxator-tk'=>1,
	EBI	    =>1,
	CLARK	    =>1,
	OneCodex    =>1,
	GOTTCHA	    =>1
	);

for (my $i=0; $i<scalar(@meth); $i++){

    if($meth[$i] =~ /(set\S{2})_(\S+)/){
	my ($n,$m)=($1,$2);

	if(defined($m) && defined($method{$m}) ){
	    $m =~ s/EBI/EBI-mg/;
	    $m =~ s/Qiime/QIIME/;

	    chomp($sen[$i]);
	    chomp($ppv[$i]);
	    printf "Lindgreen.2016\t$m\t%0.4f\t%0.4f\t%0.4f\tnotes:$n;\n", 2*($sen[$i]*$ppv[$i])/($sen[$i]+$ppv[$i]), $sen[$i], $ppv[$i];
	    
	}
	else {
	    #print "FAIL2:$meth[$i]\n";
	}

    }
    else {
	#print "FAIL1:$meth[$i]\n";
    }
}

close(IN); 




