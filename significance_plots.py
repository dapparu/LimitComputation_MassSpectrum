

import ROOT as rt
# import root_numpy as rtnp
import csv
import re
import sys
import collections
import os
sys.argv.append(' -b- ')

from collections import OrderedDict
import uproot
import pandas as pd
import math
import scipy
import awkward
import numpy as np
import time
from array import array
from histo_utilities import std_color_list, create_TGraph, find_intersect

rt.gROOT.SetBatch(True)

import CMS_lumi, tdrstyle
a = tdrstyle.setTDRStyle()
CMS_lumi.writeExtraText = 0

print(sys.version)

def setColorAndMarkerGr(gr,a,b,size=1):
    gr.SetMarkerColor(a)
    gr.SetLineColor(a)
    gr.SetMarkerStyle(b)
    gr.SetMarkerSize(size)

searchRegion='SRs'
labelSR=searchRegion+'_test'
labelSignal="ppStau"
labelSignal="gluino"
labelSignal="stop"

xsec_list = 'xSec.dat'
models = ['HSCPgluinoOnlyNeutral', 'gluino', 'gmsbStau', 'pairStau', 'stopOnlyNeutral', 'stop', 'tauPrimeCharge1e', 'tauPrimeCharge2e']
models = [ 'gluino', 'pairStau', 'stop']
if(labelSignal=="gluino"):
    models = [ 'gluino' ]
elif(labelSignal=="ppStau"):
    models = [ 'pairStau' ]
elif(labelSignal=="stop"):
    models = [ 'stop' ]


data = np.loadtxt(xsec_list, dtype=str)
xsec_hscp = {}
sample_names = []
for i in range(len(data)):
    if ('HSCP' in data[i,0]):
        name = data[i,0][:data[i,0].find('_TuneCP5')]
        xsec_hscp[name] = float(data[i, 1])
        flag = 0
        for m in models:
            if m in name:flag = 1
        if flag == 0:continue
        sample_names.append(name)

signal_names = {}
for m in models:
    mass = []
    signal_names[m] = []
    for s in sample_names:
        if m+'_' in s:
            signal_names[m].append(s)
            mass.append(int(s[s.find('M-')+2:]))
    mass = np.array(mass)
    inds = mass.argsort()
    signal_names[m] = list(np.array(signal_names[m])[inds])

print 'sample names, signal names',sample_names, signal_names


signal = []
for k,v in signal_names.items():signal += v
print 'signal',signal



limitTrees_90 =OrderedDict()
dataCards_90 = OrderedDict()
limits_90 = OrderedDict()

limitTrees_99 =OrderedDict()
dataCards_99 = OrderedDict()
limits_99 = OrderedDict()

limitTrees_999 =OrderedDict()
dataCards_999 = OrderedDict()
limits_999 = OrderedDict()

limitTrees_999muMinus2sigma =OrderedDict()
dataCards_999muMinus2sigma = OrderedDict()
limits_999muMinus2sigma = OrderedDict()


dataCardDir = '/opt/sbg/cms/safe1/cms/dapparu/HSCP/Combine/CMSSW_11_3_4/src/HiggsAnalysis/HSCPLimit/combine/datacards/mass/v_'+labelSR+'/'
limitDir = dataCardDir.replace('datacards', 'limitTrees')

dataCardDirSR1 = 'datacards_test_SR1_5may_v3/'
limitDirSR1 = 'limitTrees_test_SR1_5may_v3/'

dataCardDirSR2 = 'datacards_test_SR2_5may_v3/'
limitDirSR2 = 'limitTrees_test_SR2_5may_v3/'

dataCardDirSR3 = 'datacards_test_SR3_5may_v3/'
limitDirSR3 = 'limitTrees_test_SR3_5may_v3/'

dataCardDirSR3_muMinus2sigma = 'datacards_test_SR3_5june_muMinus2sigma/'
limitDirSR3_muMinus2sigma = 'limitTrees_test_SR3_5june_muMinus2sigma/'

for s in signal:
    limitTrees_90[s] = {}
    dataCards_90[s] = {}
    limitTrees_99[s] = {}
    dataCards_99[s] = {}
    limitTrees_999[s] = {}
    dataCards_999[s] = {}
    limitTrees_999muMinus2sigma[s] = {}
    dataCards_999muMinus2sigma[s] = {}
    if 'gluino' in s:
        name = s.replace('HSCPg', 'G')
    elif 'Stau' in s:
        name = s.replace('HSCP', '')
    elif 'stop' in s:
        name = s.replace('HSCPs', 'S')
    name = name.replace('_M-', '')
    #dataCards_90[s] = dataCardDir + '{}_nominal.txt'.format(name)
    dataCards_90[s] = dataCardDirSR1 + '{}_2018.txt'.format(name)
    dataCards_99[s] = dataCardDirSR2 + '{}_2018.txt'.format(name)
    dataCards_999[s] = dataCardDirSR3 + '{}_2018.txt'.format(name)
    dataCards_999muMinus2sigma[s] = dataCardDirSR3_muMinus2sigma + '{}_2018.txt'.format(name)
    #print dataCards_90[s]
    #print limitDir + 'higgsCombine.{}'.format(name) + '.AsymptoticLimits.mH120.root'
    limitTrees_90[s] = limitDirSR1 + 'higgsCombine.{}'.format(name) + '.Significance.mH120.root'
    limitTrees_99[s] = limitDirSR2 + 'higgsCombine.{}'.format(name) + '.Significance.mH120.root'
    limitTrees_999[s] = limitDirSR3 + 'higgsCombine.{}'.format(name) + '.Significance.mH120.root'
    limitTrees_999muMinus2sigma[s] = limitDirSR3_muMinus2sigma + 'higgsCombine.{}'.format(name) + '.Significance.mH120.root'

for i,m in enumerate(limitTrees_90.keys()):
    if not os.path.isfile(dataCards_90[m]):
        continue
    if len(uproot.open(limitTrees_90[m]).keys()) == 2:
        T = uproot.open(limitTrees_90[m])['limit']
        limits_90[m] = np.array(T.array('limit'))

for i,m in enumerate(limitTrees_99.keys()):
    if not os.path.isfile(dataCards_99[m]):
        continue
    if len(uproot.open(limitTrees_99[m]).keys()) == 2:
        T = uproot.open(limitTrees_99[m])['limit']
        limits_99[m] = np.array(T.array('limit'))

for i,m in enumerate(limitTrees_999.keys()):
    if not os.path.isfile(dataCards_999[m]):
        continue
    if len(uproot.open(limitTrees_999[m]).keys()) == 2:
        T = uproot.open(limitTrees_999[m])['limit']
        limits_999[m] = np.array(T.array('limit'))

for i,m in enumerate(limitTrees_999muMinus2sigma.keys()):
    if not os.path.isfile(dataCards_999muMinus2sigma[m]):
        continue
    if len(uproot.open(limitTrees_999muMinus2sigma[m]).keys()) == 2:
        T = uproot.open(limitTrees_999muMinus2sigma[m])['limit']
        limits_999muMinus2sigma[m] = np.array(T.array('limit'))

#print 'limits90',limits_90

# make plots



h = {}

c=rt.TCanvas()

for i, m in enumerate(signal_names.keys()):
    x = []
    ySR1 = []
    ySR2 = []
    ySR3 = []
    ySR3muMinus2sigma = []
    for key in limits_90.keys():
        if m in key and len(limits_90[key])>0:
            if (('gluino' in key) or ('Stau' in key) or ('stop' in key)):
                ySR1.append(limits_90[key][0])
            x.append(int(key[key.find('M-')+2:]))
    for key in limits_99.keys():
        if m in key and len(limits_99[key])>0:
            if (('gluino' in key) or ('Stau' in key) or ('stop' in key)):
                ySR2.append(limits_99[key][0])
            #x.append(int(key[key.find('M-')+2:]))
    for key in limits_999.keys():
        if m in key and len(limits_999[key])>0:
            if (('gluino' in key) or ('Stau' in key) or ('stop' in key)):
                ySR3.append(limits_999[key][0])
            #x.append(int(key[key.find('M-')+2:]))
    for key in limits_999muMinus2sigma.keys():
        if m in key and len(limits_999muMinus2sigma[key])>0:
            if (('gluino' in key) or ('Stau' in key) or ('stop' in key)):
                ySR3muMinus2sigma.append(limits_999muMinus2sigma[key][0])
            #x.append(int(key[key.find('M-')+2:]))

    if len(x) ==0 :continue
    #h[m+'_SR1'] = create_TGraph(x,ySR1,  axis_title=['Mass [GeV]', 'Significance'])
    #h[m+'_SR2'] = create_TGraph(x,ySR2,  axis_title=['Mass [GeV]', 'Significance'])
    h[m+'_SR3'] = create_TGraph(x,ySR3,  axis_title=['Mass [GeV]', 'Significance'])
    h[m+'_SR3muMinus2sigma'] = create_TGraph(x,ySR3muMinus2sigma,  axis_title=['Mass [GeV]', 'Significance'])

leg=rt.TLegend(0.7,0.7,0.9,0.9)

for i,m in enumerate(h.keys()):

    #h[m].SaveAs("limits_pdf_opti/significance_"+labelSignal+"_"+labelSR+".root")
    #h[m].SaveAs("limits_pdf_opti/significance_"+labelSignal+"_"+labelSR+".pdf")
    c.cd()
    if('SR3muMinus2sigma' in m):
        setColorAndMarkerGr(h[m],30,23)
        leg.AddEntry(h[m],"SR3 (#mu-2#sigma)","PE1")
    elif('SR1' in m):
        setColorAndMarkerGr(h[m],1,20)
        leg.AddEntry(h[m],"SR1","PE1")
    elif('SR2' in m):
        setColorAndMarkerGr(h[m],38,21)
        leg.AddEntry(h[m],"SR2","PE1")
    elif('SR3' in m):
        setColorAndMarkerGr(h[m],46,22)
        leg.AddEntry(h[m],"SR3","PE1")
    
    if(i==0):
        h[m].Draw("AP")
    else:
        h[m].Draw("P")
    h[m].SetMinimum(0)
    h[m].SetMaximum(15)

leg.Draw("same")
c.SetGridx()
c.SetGridy()

tdrstyle.setTDRStyle()
CMS_lumi.cmsText     = "CMS"
iPos = 0
CMS_lumi.extraText = "Internal"
CMS_lumi.writeExtraText=True
CMS_lumi.cmsText=""
CMS_lumi.extraText = "Private work"

if( iPos==0 ): CMS_lumi.relPosX = 0.12
# CMS_lumi.CMS_lumi(c, 4, 0)
CMS_lumi.lumi_13TeV  = "101 fb^{-1}"
CMS_lumi.CMS_lumi(c, 4, iPos)

c.SetTicky(1)
c.SetTickx(1)

c.Draw()

ofile_name="limits_pdf_opti/test_significance_"+labelSignal+"_"+labelSR

c.SaveAs(ofile_name+".root")
c.SaveAs(ofile_name+".pdf")
