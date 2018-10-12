#!/usr/bin/perl

use warnings;
use strict;

my $inFile = "siegwald/ranking-of-metagenomics-tools-Siegwald-cleaned.tsv";

my %method = (
    BMP      => 1,
    Kraken   => 1,
    mothur   => 1,
    OneCodex => 1,
    QIIME    => 1,
    CLARK    => 1
    );
my $lab;
    
my (%sensValues, %ppvValues, %fValues) = ((), (), ());  #hash of arrays,  

open(IN, "< $inFile");
while(my $in=<IN>){
    chomp($in);
    $in =~ s/One Codex/OneCodex/;
    $in =~ s/kraken/Kraken/;
    $in =~ s/QIIME SortMeRna SUMACLUST/QIIME/;
    $in =~ s/QIIME UCLUST/QIIME/;

#    if($in =~ /\S+\s+level;(.*error.*reads)/){
    if($in =~ /\#(\S+;\d+ reads)/){
	printMe(\%sensValues, \%ppvValues, \%fValues, $lab) if (defined($lab));
	(%sensValues, %ppvValues, %fValues) = ((), (), ()); 

	$lab=$1;
    }
    
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

	push(@{    $fValues{$aa[0]} },  $aa[4]);
	push(@{ $sensValues{$aa[0]} },  $aa[3]);
	push(@{  $ppvValues{$aa[0]} },  $aa[2]);

	push(@{    $fValues{$aa[0]} },  $aa[7]);
	push(@{ $sensValues{$aa[0]} },  $aa[6]);
	push(@{  $ppvValues{$aa[0]} },  $aa[5]);
	
	#print "Siegwald.2017\t$aa[0]\t$aa[4]\t$aa[3]\t$aa[2]\t$lab; Original values - 200(V3); Reference Database = $aa[1]; $lab2\n";
	#print "Siegwald.2017\t$aa[0]\t$aa[7]\t$aa[6]\t$aa[5]\t$lab; Original values - 400(V4-V5); Reference Database = $aa[1]; $lab2\n";
    }
}
close(IN); 

printMe(\%sensValues, \%ppvValues, \%fValues, $lab);

exit(0); 


######################################################################
#printMe:
sub printMe {
    
    my ($sensValues, $ppvValues, $fValues, $lab)=@_; 

    foreach my $tool (sort keys %{ $sensValues }){
	
	my $medSens = median( $sensValues->{$tool} ); 
	my $medPPV  = median(  $ppvValues->{$tool} ); 
	my $medF    = median(    $fValues->{$tool} ); 
	#Paper	Method	F1.measure	Sensitivity	PPV	Notes
	printf "Siegwald.2017\t$tool\t%0.3f\t%0.3f\t%0.3f\t$lab;\n", $medF, $medSens, $medPPV;
    }
    

    
}

######################################################################
#flogged from http://wiki.answers.com/Q/How_can_you_calculate_the_average_and_median_in_perl_by_subroutine
sub median { 
    @_ == 1 or die ('Sub usage: $median = median(\@array);'); 
    my ($array_ref) = @_; 
    my $count = scalar @$array_ref; 
# Sort a COPY of the array, leaving the original untouched 
    my @array = sort { $a <=> $b } @$array_ref; 
    if ($count == 1){
	return $array[0];
    }
    elsif ($count % 2) { 
	return $array[int($count/2)]; 
    } else { 
	return ($array[$count/2] + $array[$count/2 - 1]) / 2; 
    } 
} 


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

