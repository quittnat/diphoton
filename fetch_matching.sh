#!/bin/bash

REMDIR=gg_minitree_data_030903p1_14may
FILE1=matchingfile_data_14may.root
FILE2=Photon-Run2011AB-21Jun2013-v1-AOD-EXTRA-14may.root

T3DCAPHOME=dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/peruzzi/

mkdir -p /scratch/peruzzi
cd /scratch/peruzzi || exit 1
rm ${FILE1} ${FILE2}
dccp ${T3DCAPHOME}/${REMDIR}/${FILE1} .
dccp ${T3DCAPHOME}/${REMDIR}/${FILE2} .

