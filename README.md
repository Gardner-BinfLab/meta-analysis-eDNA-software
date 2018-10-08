# meta-analysis-eDNA-software

A meta-analysis of benchmarks of environmental DNA sequence analysis tools. 

##Parsing data from the following meta-genomics analysis benchmarks:

- Almeida.2018
- Bazinet.2012
- Lindgreen.2016
- McIntyre.2017
- Peabody.2015
- Sczyrba.2017
- Siegwald.2017

```
cd data/
../bin/almeida2tab.pl     >  metagenome-meta-analysis-F-measures.tsv
../bin/bazinet2tab.pl     >> metagenome-meta-analysis-F-measures.tsv
../bin/lindgreen2tab2.pl  >> metagenome-meta-analysis-F-measures.tsv
../bin/siegwald2tab2.pl   >> metagenome-meta-analysis-F-measures.tsv
```

```
##These datasets did not meet the inclusion criteria for the meta-analysis. The parsing scripts have been included for the sake of completeness:
##./bin/peabody2tab2.pl    >> data/metagenome-meta-analysis-F-measures.tsv
##cd data/sczyrba/
##paste cami-challenge.recall-precision-genus cami-challenge.recall-precision-phylum cami-challenge.method_names cami-challenge.ids | sort -k5d | perl -lane 's/Common Kmers/commonkmers/;  @F=split(/\t/, $_); printf "Sczyrba.2017\t$F[4]\t%0.4f\t%0.4f\t%0.4f\tnotes:genus:$F[5]\nSczyrba.2017\t$F[4]\t%0.4f\t%0.4f\t%0.4f\tnotes:phylum:$F[5]\n", 2*($F[0]*$F[1])/($F[0]+$F[1]), $F[0], $F[1], 2*($F[2]*$F[3])/($F[2]+$F[3]), $F[2], $F[3]'                                        >> ../metagenome-meta-analysis-F-measures.tsv
##cd -
##./bin/sczyrba2tab2.pl    >> data/metagenome-meta-analysis-F-measures.tsv
##./bin/mcintyre2tab2.pl   >> data/metagenome-meta-analysis-F-measures.tsv
```

##Plotting results:

```
./bin/visualiseResults.R
```
