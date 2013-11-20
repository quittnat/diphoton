#!/bin/bash

mkdir -p /scratch/peruzzi
cd /scratch/peruzzi || exit 1
myhost=$(hostname)
myport=500${myhost:4}
echo $myhost $myport
ssh $1 "nc $myhost $myport < $2" &
nc -l ${myport} > /scratch/peruzzi/$(basename $2)
wait