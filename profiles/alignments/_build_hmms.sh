#!/bin/sh

for i in $(ls *FASTA | sed 's/\.FASTA//g')
do 
  hmmbuild ../hmm/${i}.hmm ${i}.FASTA;
done