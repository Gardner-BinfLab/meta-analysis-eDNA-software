#!/usr/bin/perl

use warnings;
use strict;

my %method = (
    MAPseq	    =>1,
    mothur	    =>1,
    QIIME	    =>1,
    QIIME2	    =>1
    );

my $lab;

my $inFile = "almeida/Tables-combined.tsv";
open(IN, "< $inFile");


print "Paper\tMethod\tF1.measure\tSensitivity\tPPV\tNotes\n";

my (%sensValues, %ppvValues, %fValues);  #hash of arrays,  
while(my $in=<IN>){
    chomp($in);

    if($in =~ /Table S\d\.\s+(\S+)/){
	printMe(\%sensValues, \%ppvValues, \%fValues, $lab) if (defined($lab));
	(%sensValues, %ppvValues, %fValues) = ((), (), ()); 
	$lab="Dataset:$1";
	next;
    }
    
    my @a=split(/\s+/, $in);
    if(defined($a[0]) && defined($method{$a[0]}) ){
	my @aa = split(/\t/, $in);
	my $lab2 = "DB:$a[1]";
	
	#Sensitivity = TPR, Recall
	#PPV         = Precision
	#PPV = F*Sens/(2xSens - F)

	#database, taxonomy level are point estimates of the same, non-independent value
	#taxLevel: family
	my $sens = $aa[2]/100;
	my $F    = $aa[4]/100;
	my $ppv  = $F*$sens / (2*$sens - $F); #PPV values not provided, so calculated them...
	push(@{    $fValues{$aa[0]} },  $F);
	push(@{ $sensValues{$aa[0]} },  $sens);
	push(@{  $ppvValues{$aa[0]} },  $ppv);

	#taxLevel: genus
	$sens = $aa[6]/100;
	$F    = $aa[8]/100;
	$ppv  = ($F*$sens / (2*$sens - $F));
	push(@{    $fValues{$aa[0]} },  $F);
	push(@{ $sensValues{$aa[0]} },  $sens);
	push(@{  $ppvValues{$aa[0]} },  $ppv);
	
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
	printf "Almeida.2018\t$tool\t%0.3f\t%0.3f\t%0.3f\t$lab;\n", $medF, $medSens, $medPPV;
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
