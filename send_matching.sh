for i in `seq 30 39`; do qsub -q all.q -l hostname=t3wn$i fetch_matching.sh; done
