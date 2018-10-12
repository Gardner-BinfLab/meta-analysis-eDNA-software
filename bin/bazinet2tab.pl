#!/usr/bin/perl

use warnings;
use strict;

my $inFile = 'bazinet/ranking-of-metagenomics-tools-Bazinet.tsv';

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

#print "Paper\tMethod\tF1.measure\tSensitivity\tPPV\tNotes\n";
open(IN, "< $inFile");
my $lab;
my (%sensValues, %ppvValues, %fValues) = ((), (), ());  #hash of arrays,  
while(my $in=<IN>){
    chomp($in);

    $in =~ s/CARMA/CARMA3/g;
    $in =~ s/PhylopythiaS/PhyloPythiaS/g;
    $in =~ s/PhymmBL/phymmBL/g;
    $in =~ s/Ralphy/RAIphy/g;
    $in =~ s/RAlphy/RAIphy/g;
    
    if($in =~ /(.*dataset)/){
	#printMe(\%sensValues, \%ppvValues, \%fValues, $lab) if (defined($lab));
	#(%sensValues, %ppvValues, %fValues) = ((), (), ()); 
	
	$lab=$1;
    }
    
    my @a=split(/\s+/, $in);
    my @aa = split(/\t/, $in);
    #print "a:[[[@a]]]\n";
    if(defined($a[0]) && defined($method{$a[0]}) && defined($aa[3]) && defined($aa[2]) && isNumeric($aa[3]) && isNumeric($aa[2]) ){
	my($sens, $ppv) = ($aa[2]/100, $aa[3]/100);
	my $f1          = (2*$sens*$ppv)/($sens+$ppv);
	printf "Bazinet.2012\t$aa[0]\t%0.4f\t%0.4f\t%0.4f\t$lab; \n", $f1, $sens, $ppv;
    }
}
close(IN); 
#printMe(\%sensValues, \%ppvValues, \%fValues, $lab);


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
	printf "Bazinet.2012\t$tool\t%0.3f\t%0.3f\t%0.3f\t$lab;\n", $medF, $medSens, $medPPV;
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
