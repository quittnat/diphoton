#!/usr/bin/python

##################################################################

file=[]
modes=[]
isdata=[]

# DATA (do not swap with mc)
file.append('./gg_minitree_data_030903p1_06apr/Photon-Run2011AB-21Jun2013-v1-AOD.root')
modes.append(['Default','ESCALEup','ESCALEdown','JECup','JECdown'])
isdata.append(1)

# MC (do not swap with data)
file.append('./gg_minitree_mc_030903p1_06apr/sig.root')
modes.append(['Default','ESMEARup','ESMEARdown','JERup','JERdown'])
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


for i in xrange(2):
    if (isdata[i]==0):
        strdata='sig'
    else:
        strdata='data'
    for mode in modes[i]:
        wait_processes()
        thisnumber = number
        if (mode.rfind('standard')>=0):
            thisnumber = -1
        args = ['root','-q','-b','-l','template_production.C+O("'+file[i]+'","effunf",'+str(isdata[i])+',"outphoton_effunf_'+strdata+'_'+mode+'.root","photoniso",'+str(thisnumber)+',false,"","'+mode+'");']
        #    print 'Running root'+args
        lista_processi.append(Popen(args))
            
for i in lista_processi:
    i.wait()


