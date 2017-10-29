#!/usr/bin/perl

use warnings;
use strict;

#my $inFile = "mcintyre2017-supp-table.txt";
my $inFile = "benchmark-papers/mcintyre-genomebiology-supplement/13059_2017_1299_MOESM9_ESM.csv";

# #benchmark-papers/mcintyre-genomebiology-supplement/13059_2017_1299_MOESM9_ESM.csv
# head -n 1 benchmark-papers/mcintyre-genomebiology-supplement/13059_2017_1299_MOESM9_ESM.csv | tr "\t" "\n" | nl
#      1	Dataset
#      2	Metric
#      3	CLARK
#      4	CLARK-S
#      5	Kraken
#      6	LMAT
#      7	BlastMegan
#      8	DiamondMegan
#      9	NBC
#     10	Meta-Classifier (Precision)
#     11	Meta-Classifier (Recall)

# Dataset	Metric
# HC1	Precision
# 	Recall
# HC2	Precision
# 	Recall

my %inclMethods = (
    CLARK      => 2,
    'CLARK-S'  => 3,
    Kraken     => 4,
    LMAT       => 5,
    MEGAN      => 6,
    'MEGAN-D'  => 7,
    NBC        => 8
);

my (%precisionMethods, %recallMethods);

open(IN, "< $inFile");
my $dataset = '';
my $lab='';
my ($methCnt, $datasetCnt)=(0,0);
while(my $in=<IN>){
    next if ($in =~ /^Dataset/);
    chomp($in);
    
    my @a=split(/\s+/, $in);
#    $in =~ s/\t/|/g; 
#    print "$in\n"; 
    
    if($in =~ /(\S+)\s+Precision/){
	$dataset = $1;
	foreach my $m (keys %inclMethods){
	    $precisionMethods{ $m } = $a[$inclMethods{$m}] if defined( $a[$inclMethods{$m}] );
	}
    }
    elsif($in =~ /Recall/){
	foreach my $m (sort keys %inclMethods){
	    $recallMethods{ $m } = $a[$inclMethods{$m}] if defined( $a[$inclMethods{$m}] );
	    
	    my ($sens,$ppv) = ($recallMethods{$m}/100, $precisionMethods{$m}/100);
	    
	    printf "McIntyre.2017\t$m\t%0.4f\t%0.4f\t%0.4f\tTableS7:$dataset;Precision:$precisionMethods{$m};Recall:$recallMethods{$m};\n", 2*($sens*$ppv)/($sens + $ppv), $sens, $ppv if (defined($precisionMethods{$m}) and defined($recallMethods{$m}));
	}
	$datasetCnt++;
	$dataset = "";
    }
    
}

close(IN); 


######################################################################

$inFile = "benchmark-papers/mcintyre-genomebiology-supplement/13059_2017_1299_MOESM4_ESM.csv";

# #benchmark-papers/mcintyre-genomebiology-supplement/13059_2017_1299_MOESM4_ESM.csv    

# head -n 2 benchmark-papers/mcintyre-genomebiology-supplement/13059_2017_1299_MOESM4_ESM.csv | tail -n 1 | tr "\t" "\n" | nl
#      1	dataset
#      2	measure
#      3	BlastMegan_filtered
#      4	BlastMegan_filtered_liberal
#      5	CLARK
#      6	CLARK-S
#      7	DiamondMegan_filtered
#      8	GOTTCHA
#      9	Kraken
#     10	Kraken_filtered
#     11	LMAT
#     12	MetaFlow
#     13	MetaPhlAn
#     14	NBC
#     15	PhyloSift
#     16	PhyloSift_filtered

# dataset	measure
# BioPool	precision
# BioPool	recall
# HC1	precision
# HC1	recall
# HC2	precision
# HC2	recall

#only using unfiltered results (i.e. defaults):
%inclMethods = (
    MEGAN          =>  2,
    CLARK          =>  4,
    'CLARK-S'      =>  5,
    'MEGAN-D'      =>  6,
    GOTTCHA        =>  7,
    Kraken         =>  8,
    LMAT           => 10,
    MetaFlow       => 11,
    'MetaPhlAn2.0' => 12,
    NBC            => 13,
    PhyloSift      => 14
);

(%precisionMethods, %recallMethods) = ((),());

open(IN, "< $inFile");
while(my $in=<IN>){
    #print $in; 
    next if ($in =~ /^dataset/ or $in =~ /^Table/);
    chomp($in);
    
    my @a=split(/\s+/, $in);
    if($in =~ /(\S+)\s+precision/){
	$dataset = $1;
	foreach my $m (keys %inclMethods){
	    $precisionMethods{ $m } = $a[$inclMethods{$m}] if defined( $a[$inclMethods{$m}] );
	}
    }
    elsif($in =~ /recall/){
	foreach my $m (sort keys %inclMethods){
	    $recallMethods{ $m } = $a[$inclMethods{$m}] if defined( $a[$inclMethods{$m}] );
	    my ($sens,$ppv) = ($recallMethods{$m}/100, $precisionMethods{$m}/100);
	    my $sum = $sens+$ppv;
	    $sum = 1.0 if ($sum==0.0);
	    printf "McIntyre.2017\t$m\t%0.4f\t%0.4f\t%0.4f\tTableS3:$dataset;Precision:$precisionMethods{$m};Recall:$recallMethods{$m};\n", 2*($sens*$ppv)/($sum), $sens, $ppv if (defined($precisionMethods{$m}) and defined($recallMethods{$m}));
	}
	$datasetCnt++;
	$dataset = "";
    }
    
}

close(IN); 



