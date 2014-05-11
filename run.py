#!/usr/bin/python

##################################################################

file='./gg_minitree_data_030903p1_06apr/Photon-Run2011AB-21Jun2013-v1-AOD.root'
modes=['sigsig','sigbkg','bkgbkg','standard']
isdata=1
number=5e5


#file='./gg_minitree_mc_030903p1_22mar/allmc.root'
#modes=['signal','randomcone','background','sieiesideband','sigsig','sigbkg','bkgbkg','standard']
#isdata=0
#number=5e5

##################################################################

spawn=4
lista_processi=[]

from subprocess import Popen
from time import sleep

def wait_processes():
    sleep(1)
    for i in lista_processi:
        if (i.poll()!=None):
            lista_processi.remove(i)
    if (len(lista_processi)>=spawn):
        wait_processes()
    else:
        return

if (isdata==0):
    strdata='allmc'
else:
    strdata='data'

for mode in modes:
    wait_processes()
    thisnumber = number
    if (mode.rfind('standard')>=0):
        thisnumber = -1
    args = ['root','-q','-b','-l','template_production.C+O("'+file+'","'+mode+'",'+str(isdata)+',"outphoton_'+strdata+'_'+mode+'.root","photoniso",'+str(thisnumber)+');']
#    print 'Running root'+args
    lista_processi.append(Popen(args))

for i in lista_processi:
    i.wait()




