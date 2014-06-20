#!/usr/bin/python

##################################################################

file=[]
modes=[]
isdata=[]

# DO NOT CHANGE THE ORDER OF WHAT FOLLOWS

# DATA (do not swap with mc)
file.append('./gg_minitree_data_030903p1_17jun14/Photon-Run2011AB-21Jun2013-v1-AOD.root')
modes.append(['Default','ESCALEup','ESCALEdown','JECup','JECdown'])
isdata.append(1)

# MC (do not swap with data)
file.append('./gg_minitree_mc_030903p1_17jun14_lightforPU/sig.root')
modes.append(['Default','ESMEARup','ESMEARdown','JERup','JERdown'])
isdata.append(0)

# MC pileup rew up/down
file.append('./gg_minitree_mc_030903p1_17jun14_lightforPU_UPvar/sig.root')
modes.append(['Default'])
isdata.append(0)
file.append('./gg_minitree_mc_030903p1_17jun14_lightforPU_DOWNvar/sig.root')
modes.append(['Default'])
isdata.append(0)


#################################################################

number=-1
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


for i in xrange(4):
    if (isdata[i]==0):
        strdata='sig'
    else:
        strdata='data'
    for mode in modes[i]:
        wait_processes()
        thisnumber = number
        tag=mode
        if (i==2):
            tag='PUup'
        if (i==3):
            tag='PUdown'
        args = ['root','-q','-b','-l','template_production.C+O("'+file[i]+'","effunf",'+str(isdata[i])+',"outphoton_effunf_'+strdata+'_'+tag+'.root","photoniso",'+str(thisnumber)+',false,"","'+mode+'");']
        #    print 'Running root'+args
        lista_processi.append(Popen(args))
            
for i in lista_processi:
    i.wait()


