#Set MAX to listauuid.size()/20 approx. by defect
MAX=143
for i in $(seq 0 ${MAX}); do qsub -q all.q -l hostname=t3wn3[0-9] run_prepare_matchingfile_forstep2.sh $i; done
