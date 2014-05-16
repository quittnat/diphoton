#!/usr/bin/python
from ROOT import *

file=TFile('outphoton_allmc_sigbkg.root')
file.cd('roofit')

out=[]

for reg in ['EBEB','EBEE','EEEE']:
    t1=gDirectory.Get('template_roodset_'+reg+'_bkgsig')
    t2=gDirectory.Get('template_roodset_'+reg+'_sigbkg')
    t1.Draw("rooissigcont>>htemp","rooweight")
    t2.Draw("rooissigcont>>htemp2","rooweight")
    if (reg=='EBEE'):
        out.append(htemp.GetMean())
        out.append(htemp2.GetMean())
    else:
        out.append((htemp.GetMean()*htemp.Integral()+htemp2.GetMean()*htemp2.Integral())/(htemp.Integral()+htemp2.Integral()))

print '### Add this in binsdef.h'    
print 'const double sig_contamination_in_templates[4] = {'+','.join([str(x) for x in out])+'};'
