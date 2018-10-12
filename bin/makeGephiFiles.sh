#!/bin/sh

echo "from,to,type,weight" > medianVals-gephiEdges.tsv
cat medianVals.tsv | sort -k1,1d -k3,3nr | perl -lane 'if(defined($prevM) && ($prevP eq $F[0]) ){$cnt++; print "$F[1],$prevM,$F[0],$cnt"}else{$cnt=0;} ($prevP, $prevM)=($F[0], $F[1]); ' >> medianVals-gephiEdges.tsv



