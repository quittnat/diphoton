for i in `seq 30 39`; do qsub -q short.q -l hostname=t3wn$i fetch_matching.sh $1 $2; done
