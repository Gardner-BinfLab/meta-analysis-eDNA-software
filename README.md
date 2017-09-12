# meta-analysis-eDNA-software

A meta-analysis of benchmarks of environmental DNA sequence analysis tools. 


######################################################################
#Parsing data from the following manuscripts:
# Bazinet.2012
# Lindgreen.2016
# Siegwald.2017
# Peabody.2015
# Sczyrba.2017
# McIntyre.2017

./bin/bazinet2tab.pl      > data/metagenome-meta-analysis-F-measures.tsv
./bin/lindgreen2tab2.pl  >> data/metagenome-meta-analysis-F-measures.tsv
./bin/siegwald2tab2.pl   >> data/metagenome-meta-analysis-F-measures.tsv
./bin/peabody2tab2.pl    >> data/metagenome-meta-analysis-F-measures.tsv
cd data/sczyrba/
paste cami-challenge.recall-precision-genus cami-challenge.recall-precision-phylum cami-challenge.method_names cami-challenge.ids | sort -k5d | perl -lane 's/Common Kmers/commonkmers/;  @F=split(/\t/, $_); printf "Sczyrba.2017\t$F[4]\t%0.4f\t%0.4f\t%0.4f\tnotes:genus:$F[5]\nSczyrba.2017\t$F[4]\t%0.4f\t%0.4f\t%0.4f\tnotes:phylum:$F[5]\n", 2*($F[0]*$F[1])/($F[0]+$F[1]), $F[0], $F[1], 2*($F[2]*$F[3])/($F[2]+$F[3]), $F[2], $F[3]'                                        >> ../metagenome-meta-analysis-F-measures.tsv
cd -
./bin/sczyrba2tab2.pl    >> data/metagenome-meta-analysis-F-measures.tsv
./bin/mcintyre2tab2.pl   >> data/metagenome-meta-analysis-F-measures.tsv



######################################################################
#Plotting results:

./bin/visualiseResults.R

