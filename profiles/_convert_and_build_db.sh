#/bin/sh

for i in $(ls *hmm | sed 's/\.hmm//g')
do 
  hmmconvert ${i}.hmm > ${i}_new.hmm;
done

cat *_new.hmm >> database.hmm

rm *_new.hmm