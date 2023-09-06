"""

2022/10/03
create datacard from the mass distributions of signal and data in region D
using distributions that are already normalized from Dylan

"""

import ROOT as rt
import csv
import re
import sys
import collections
import os
sys.argv.append(' -b- ')

from collections import OrderedDict
import uproot

import scipy
import awkward
import numpy as np
import time

from histo_utilities import std_color_list, create_TGraph, find_intersect


rt.gROOT.SetBatch(True)

import CMS_lumi, tdrstyle
a = tdrstyle.setTDRStyle()
CMS_lumi.writeExtraText = 0

def setColorAndMarkerGr(gr,a,b):
    gr.SetMarkerColor(a)
    gr.SetLineColor(a)
    gr.SetMarkerStyle(b)
    gr.SetMarkerSize(2)

def statErr(h1,name):
    statErr=h1.Clone()
    statErr.Reset()
    statErr.SetName(name)
    for i in range (1,statErr.GetNbinsX()):
        if h1.GetBinContent(i)>0:
            #print (i,statErr.GetBinError(i),statErr.GetBinContent(i),statErr.GetBinError(i)/statErr.GetBinContent(i))
            statErr.SetBinContent(i,h1.GetBinError(i)+h1.GetBinContent(i))
            #print i,h1.GetBinLowEdge(i), h1.GetBinError(i),h1.GetBinContent(i), h1.GetBinError(i)/h1.GetBinContent(i)
            #statErr.SetBinContent(i,statErr.GetBinError(i))
        else:
            statErr.SetBinContent(i,1)
    return statErr

def BiasCorrection(h1,a_,b_):
    h=h1.Clone()
    for i in range (0,h.GetNbinsX()+1):
        mass = h.GetBinLowEdge(i)
        if(mass<25): continue
        h.SetBinContent(i,h.GetBinContent(i)*(a_*mass+b_))
    return h

def make_datacard_hscp(outDataCardsDir,  modelName, signal, bkg, observation, sig_unc, bkg_unc):


    text_file = open(outDataCardsDir+modelName+".txt", "w")

    text_file.write('imax {0} \n'.format(1))
    text_file.write('jmax {0} \n'.format(1))
    text_file.write('kmax * \n')
    text_file.write('shapes * * FAKE \n')



    text_file.write('--------------- \n')
    text_file.write('--------------- \n')
    text_file.write('bin \t  chD \n')
    text_file.write('observation \t {0:6.2f} \n'.format(observation))
    text_file.write('------------------------------ \n')
    text_file.write('bin \t chD \t chD \n')
    text_file.write('process \t signal \t bkg \n')
    text_file.write('process \t 0 \t 1 \t \n')
    text_file.write('rate \t {} \t {} \n'.format(signal, bkg))
    text_file.write('------------------------------ \n')

    bckg_up=0
    bckg_down=0
    signal_up=0
    signal_down=0
#   #### uncertainties ####
    for k,v in sig_unc.items():
        if len(v)==2:
            text_file.write('{} \t lnN \t {}/{} \t - \n'.format(k, v[0],v[1]))
            signal_down=v[0]
            signal_up=v[1]
        else:text_file.write('{} \t lnN \t {} \t - \n'.format(k, v[0]))
        #text_file.write('{} \t lnN \t {} \t - \n'.format(k, 1+v))

    for k,v in bkg_unc.items():
        if len(v)==2:
            text_file.write('{} \t lnN \t -  \t {}/{} \n'.format(k, v[0],v[1]))
            bckg_down=v[0]
            bckg_up=v[1]
        else:text_file.write('{} \t lnN \t -  \t {} \n'.format(k, v[0]))

    text_file.close()

def make_datacard_hscp_combining2017and2018(outDataCardsDir,  modelName, signal2017, signal2018, bkg2017, bkg2018, observation2017, observation2018, sig_2017_unc, sig_2018_unc, sig_unc_correlated, bkg_2017_unc, bkg_2018_unc, bkg_unc_correlated):

    text_file = open(outDataCardsDir+modelName+".txt", "w")

    text_file.write('imax {0} \n'.format(2))
    text_file.write('jmax {0} \n'.format(1))
    text_file.write('kmax * \n')
    text_file.write('shapes * * FAKE \n')

    text_file.write('--------------- \n')
    text_file.write('--------------- \n')
    text_file.write('bin \t  Ch2017 \t Ch2018 \n')
    text_file.write('observation \t {} \t {} \n'.format(observation2017,observation2018))
    text_file.write('------------------------------ \n')
    text_file.write('bin \t Ch2017 \t Ch2017 \t Ch2018 \t Ch2018 \n')
    text_file.write('process \t signal \t bkg \t signal \t bkg \n')
    text_file.write('process \t 0 \t 1 \t 0 \t 1 \t \n')
    text_file.write('rate \t {} \t {} \t {} \t {} \n'.format(signal2017, bkg2017, signal2018, bkg2018))
    text_file.write('------------------------------ \n')

#   #### uncertainties ####
    for k,v in sig_2017_unc.items():
        if len(v)==2:
            text_file.write('{} \t lnN \t {}/{} \t - \t - \t - \n'.format(k, v[0],v[1]))
        else:text_file.write('{} \t lnN \t {} \t - \t - \t - \n'.format(k, v[0]))

    for k,v in sig_2018_unc.items():
        if len(v)==2:
            text_file.write('{} \t lnN \t - \t - \t {}/{} \t - \n'.format(k, v[0],v[1]))
        else:text_file.write('{} \t lnN \t - \t - \t {} \t - \n'.format(k, v[0]))

    for k,v in sig_unc_correlated.items():
        if len(v)==4:
            text_file.write('{} \t lnN \t {}/{} \t - \t {}/{} \t - \n'.format(k, v[0],v[1], v[2],v[3]))
        else:text_file.write('{} \t lnN \t {} \t - \t {} \t - \n'.format(k, v[0], v[1]))

    for k,v in bkg_2017_unc.items():
        if len(v)==2:
            vPrime=max(abs(1-v[0]),abs(1-v[1]))+1
            text_file.write('{} \t lnN \t -  \t {} \t - \t - \n'.format(k, vPrime))
        else:text_file.write('{} \t lnN \t -  \t {} \t - \t - \n'.format(k, v[0]))

    for k,v in bkg_2018_unc.items():
        if len(v)==2:
            vPrime=max(abs(1-v[0]),abs(1-v[1]))+1
            text_file.write('{} \t lnN \t - \t - \t - \t {} \n'.format(k, vPrime))
        else:text_file.write('{} \t lnN \t - \t - \t - \t {} \n'.format(k, v[0]))

    for k,v in bkg_unc_correlated.items():
        if len(v)==4:
            vPrime=max(abs(1-v[0]),abs(1-v[1]))+1
            vSecond=max(abs(1-v[2]),abs(1-v[3]))+1
            text_file.write('{} \t lnN \t - \t {} \t - \t {} \n'.format(k, vPrime, vSecond))
        else:text_file.write('{} \t lnN \t - \t {} \t - \t {} \n'.format(k, v[0], v[1]))

    #for k,v in bkg_2017_unc.items():
    #    if len(v)==2:
    #        text_file.write('{} \t lnN \t -  \t {}/{} \t - \t - \n'.format(k, v[0],v[1]))
    #    else:text_file.write('{} \t lnN \t -  \t {} \t - \t - \n'.format(k, v[0]))
#
    #for k,v in bkg_2018_unc.items():
    #    if len(v)==2:
    #        text_file.write('{} \t lnN \t - \t - \t - \t {}/{} \n'.format(k, v[0],v[1]))
    #    else:text_file.write('{} \t lnN \t - \t - \t - \t {} \n'.format(k, v[0]))
#
    #for k,v in bkg_unc_correlated.items():
    #    if len(v)==4:
    #        text_file.write('{} \t lnN \t - \t {}/{} \t - \t {}/{} \n'.format(k, v[0],v[1], v[2],v[3]))
    #    else:text_file.write('{} \t lnN \t - \t {} \t - \t {} \n'.format(k, v[0], v[1]))

    text_file.close()

def totalUncertainy(unc,i=1):
    total=0
    for k,v in unc.items():
        if len(v)==2 and (i==0 or i==1):
            total+=pow(1-v[i],2)
        else:total+=pow(1-v[0],2)
    return np.sqrt(total)

def makeYieldFile(text_file_tex,modelName,bkg,bckg_up,bckg_down,signal,signal_up,signal_down):
    #typeOfDisplay=""
    #print bkg,bckg_up,bckg_down,bkg*bckg_up,bkg*bckg_down
    typeOfDisplay='.2E'
    text_file_tex.write('\n '+modelName+' & $'+str(format(bkg,typeOfDisplay))+'^{+'+str(format(bkg*bckg_up,typeOfDisplay))+'}_{-'+str(format(bkg*bckg_down,typeOfDisplay))+'}$ & $'+str(format(signal,typeOfDisplay))+'^{+'+str(format(signal*signal_up,typeOfDisplay))+'}_{-'+str(format(signal*signal_down,typeOfDisplay))+'}$ \\\\')
    #text_file_tex.write('\n '+modelName+' & $'+str(bkg)+'^{+'+str(bckg_up)+'}_{-'+str(bckg_down)+'}$ & $'+str(signal)+'^{+'+str(signal_up)+'}_{-'+str(signal_down)+'}$ \\\\')
    text_file_tex.write('\n \hline')  
    
def fillH2(h2,targetMass,mean,stddev,s):
    #targetMass/=1000
    #mean/=1000
    #stddev/=1000
    #print targetMass, mean, stddev
    xmin=mean-stddev
    if (xmin<300): 
        xmin=300
    xmax=100000
    for i in range(0,h2.GetNbinsX()+1):
        if (i!=h2.GetXaxis().FindBin(targetMass)): 
                continue
        else:
            #for j in range(h2.GetYaxis().FindBin(mean-stddev),h2.GetYaxis().FindBin(mean+2*stddev)):
            for j in range(h2.GetYaxis().FindBin(xmin),h2.GetYaxis().FindBin(xmax)):
                if("Gluino" in s):
                    h2.SetBinContent(i,j,1)
                elif("pairStau" in s):
                    h2.SetBinContent(i,j,2)
                elif("Stop" in s):
                    h2.SetBinContent(i,j,3)

def pushSyst(syst1,syst2):
    syst1=abs(1-syst1)
    syst2=abs(1-syst2)
    res=max(syst1,syst2)
    return res*100

def integralHisto(h,xmin,xmax):
    return h.Integral(h.FindBin(xmin),h.FindBin(xmax))



if __name__ == '__main__':

    f=rt.TF1("f","exp(-x/300)",0,1000)
    integralF = f.Integral(30,1000)/f.Integral(0,1000)
    #print integralF

# # load mass distributions
    fpath =OrderedDict()
    fpathPred =OrderedDict()
    massPlotsSignal =OrderedDict()
    mass = OrderedDict()
    mass_plot = OrderedDict()

    regionSignal='SR3'
    regionBckg=''
    systSignal="Stau"
    
    if(regionSignal=='SR1'):
        regionBckg='90ias100'
    if(regionSignal=='SR2'):
        regionBckg='99ias100'
    if(regionSignal=='SR3'):
        regionBckg='999ias100'

    text_file_tex = open('yieldDir/yield'+regionSignal+'_2017_5may.tex', "w")
    text_file_tex.write('\n \\documentclass{article}')
    text_file_tex.write('\n \\begin{document}')
    text_file_tex.write('\n \\begin{center}')
    text_file_tex.write('\n \\begin{tabular}{ |l|c|c| } ')
    text_file_tex.write('\n \hline')
    text_file_tex.write('\n Yield '+regionSignal+' & Pred. & Signal \\\\')
    text_file_tex.write('\n \hline')
    text_file_tex.write('\n \hline')

    text_file_tex_2018 = open('yieldDir/yield'+regionSignal+'_2018_5may.tex', "w")
    text_file_tex_2018.write('\n \\documentclass{article}')
    text_file_tex_2018.write('\n \\begin{document}')
    text_file_tex_2018.write('\n \\begin{center}')
    text_file_tex_2018.write('\n \\begin{tabular}{ |l|c|c| } ')
    text_file_tex_2018.write('\n \hline')
    text_file_tex_2018.write('\n Yield '+regionSignal+' & Pred. & Signal \\\\')
    text_file_tex_2018.write('\n \hline')
    text_file_tex_2018.write('\n \hline')

    # load root file
    path = "/opt/sbg/cms/safe1/cms/dapparu/HSCP/CMSSW_10_6_27/src/SUSYBSMAnalysis/BackgroundPrediction/"
    pathSignal = "/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/"
    pathPred = "/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/"

    yearSignal='2018'
    codeVersionSignal='73p3_Signal'
    codeVersionSignal='77p1_Signal'

    ofileBase = rt.TFile("base.root","RECREATE")


    #fpath['Gluino400_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPgluino_M-400_CodeV'+codeVersionSignal+'_v1.root'
    fpath['Gluino500_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPgluino_M-500_CodeV'+codeVersionSignal+'_v1.root'
    fpath['Gluino800_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPgluino_M-800_CodeV'+codeVersionSignal+'_v1.root'
    fpath['Gluino1000_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPgluino_M-1000_CodeV'+codeVersionSignal+'_v1.root'
    fpath['Gluino1400_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPgluino_M-1400_CodeV'+codeVersionSignal+'_v1.root'
    fpath['Gluino1600_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPgluino_M-1600_CodeV'+codeVersionSignal+'_v1.root'
    fpath['Gluino1800_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPgluino_M-1800_CodeV'+codeVersionSignal+'_v1.root'
    fpath['Gluino2000_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPgluino_M-2000_CodeV'+codeVersionSignal+'_v1.root'
    fpath['Gluino2200_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPgluino_M-2200_CodeV'+codeVersionSignal+'_v1.root'
    fpath['Gluino2400_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPgluino_M-2400_CodeV'+codeVersionSignal+'_v1.root'
    fpath['Gluino2600_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPgluino_M-2600_CodeV'+codeVersionSignal+'_v1.root'
    fpath['Stop800_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPstop_M-800_CodeV'+codeVersionSignal+'_v1.root'
    fpath['Stop1000_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPstop_M-1000_CodeV'+codeVersionSignal+'_v1.root'
    fpath['Stop1200_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPstop_M-1200_CodeV'+codeVersionSignal+'_v1.root'
    fpath['Stop1400_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPstop_M-1400_CodeV'+codeVersionSignal+'_v1.root'
    fpath['Stop1600_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPstop_M-1600_CodeV'+codeVersionSignal+'_v1.root'
    fpath['Stop1800_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPstop_M-1800_CodeV'+codeVersionSignal+'_v1.root'
    fpath['Stop2000_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPstop_M-2000_CodeV'+codeVersionSignal+'_v1.root'
    fpath['Stop2200_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPstop_M-2200_CodeV'+codeVersionSignal+'_v1.root'
    fpath['Stop2400_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPstop_M-2400_CodeV'+codeVersionSignal+'_v1.root'
    fpath['Stop2600_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPstop_M-2600_CodeV'+codeVersionSignal+'_v1.root'
    #fpath['pairStau247_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPpairStau_M-247_CodeV'+codeVersionSignal+'_v1.root'
    fpath['pairStau308_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPpairStau_M-308_CodeV'+codeVersionSignal+'_v1.root'
    fpath['pairStau432_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPpairStau_M-432_CodeV'+codeVersionSignal+'_v1.root'
    fpath['pairStau557_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPpairStau_M-557_CodeV'+codeVersionSignal+'_v1.root'
    fpath['pairStau651_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPpairStau_M-651_CodeV'+codeVersionSignal+'_v1.root'
    fpath['pairStau745_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPpairStau_M-745_CodeV'+codeVersionSignal+'_v1.root'
    fpath['pairStau871_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPpairStau_M-871_CodeV'+codeVersionSignal+'_v1.root'
    fpath['pairStau1029_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPpairStau_M-1029_CodeV'+codeVersionSignal+'_v1.root'
    #fpath['DYcharge1e_100'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge1e_M-100_CodeV80p0_v1.root'
    #fpath['DYcharge1e_200'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge1e_M-200_CodeV80p0_v1.root'
    #fpath['DYcharge1e_400'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge1e_M-400_CodeV80p0_v1.root'
    fpath['DYcharge1e_500_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge1e_M-500_CodeV80p0_v1.root'
    fpath['DYcharge1e_800_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge1e_M-800_CodeV80p0_v1.root'
    fpath['DYcharge1e_1000_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge1e_M-1000_CodeV80p0_v1.root'
    fpath['DYcharge1e_1400_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge1e_M-1400_CodeV80p0_v1.root'
    fpath['DYcharge1e_1800_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge1e_M-1800_CodeV80p0_v1.root'
    fpath['DYcharge1e_2200_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge1e_M-2200_CodeV80p0_v1.root'
    fpath['DYcharge1e_2600_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge1e_M-2600_CodeV80p0_v1.root'
    #fpath['DYcharge2e_100'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge2e_M-100_CodeV80p0_v1.root'
    #fpath['DYcharge2e_200'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge2e_M-200_CodeV80p0_v1.root'
    fpath['DYcharge2e_400_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge2e_M-400_CodeV80p0_v1.root'
    fpath['DYcharge2e_500_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge2e_M-500_CodeV80p0_v1.root'
    fpath['DYcharge2e_800_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge2e_M-800_CodeV80p0_v1.root'
    fpath['DYcharge2e_1000_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge2e_M-1000_CodeV80p0_v1.root'
    fpath['DYcharge2e_1400_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge2e_M-1400_CodeV80p0_v1.root'
    fpath['DYcharge2e_1800_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge2e_M-1800_CodeV80p0_v1.root'
    fpath['DYcharge2e_2200_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge2e_M-2200_CodeV80p0_v1.root'
    fpath['DYcharge2e_2600_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge2e_M-2600_CodeV80p0_v1.root'

    idirSignal='HSCParticleAnalyzer/BaseName/'
    #searchRegion='SR2'
    searchRegion=regionSignal
    regionBckg=''

    if(searchRegion=='SR1'):
        regionBckg='90ias100'
    if(searchRegion=='SR2'):
        regionBckg='99ias100'
    if(searchRegion=='SR3'):
        regionBckg='999ias100'

    massPlotsSignal['Nominal'] = idirSignal+'PostS_'+searchRegion+'_Mass'
    massPlotsSignal['PU_up'] = idirSignal+'PostS_'+searchRegion+'_Mass_Pileup_up'
    massPlotsSignal['PU_down'] = idirSignal+'PostS_'+searchRegion+'_Mass_Pileup_down'
    massPlotsSignal['Fpix_up'] = idirSignal+'PostS_'+searchRegion+'_Mass_ProbQNoL1_up'
    massPlotsSignal['Fpix_down'] = idirSignal+'PostS_'+searchRegion+'_Mass_ProbQNoL1_down'
    massPlotsSignal['Gstrip_up'] = idirSignal+'PostS_'+searchRegion+'_Mass_Ias_up'
    massPlotsSignal['Gstrip_down'] = idirSignal+'PostS_'+searchRegion+'_Mass_Ias_down'
    massPlotsSignal['Pt_up'] = idirSignal+'PostS_'+searchRegion+'_Mass_Pt_up'
    massPlotsSignal['Pt_down'] = idirSignal+'PostS_'+searchRegion+'_Mass_Pt_down'
    massPlotsSignal['Trigger_up'] = idirSignal+'PostS_'+searchRegion+'_Mass_Trigger_up'
    massPlotsSignal['Trigger_down'] = idirSignal+'PostS_'+searchRegion+'_Mass_Trigger_down'
    massPlotsSignal['K_up'] = idirSignal+'PostS_'+searchRegion+'_Mass_K_up1'
    massPlotsSignal['K_down'] = idirSignal+'PostS_'+searchRegion+'_Mass_K_down1'
    massPlotsSignal['C_up'] = idirSignal+'PostS_'+searchRegion+'_Mass_C_up1'
    massPlotsSignal['C_down'] = idirSignal+'PostS_'+searchRegion+'_Mass_C_down1'

    year='2017'
    codeVersion='73p3_v4'
    nPE='200'
    endLabel='_19april'
    if(searchRegion=='SR3'):
        #endLabel='_21april_SR3'
        endLabel='_21aug'

    codeVersion='UnB_v1_v1'
    endLabel='_UnB_v3_Data_v1'

    #fpathPred['pred_2017_nominal'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    #fpathPred['pred_2017_etaup'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta2_rebinIh4_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    #fpathPred['pred_2017_etadown'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta8_rebinIh4_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    #fpathPred['pred_2017_ihup'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh2_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    #fpathPred['pred_2017_ihdown'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh8_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    #fpathPred['pred_2017_momup'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP1_rebinMass1_nPE'+nPE+endLabel+'.root')
    #fpathPred['pred_2017_momdown'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP4_rebinMass1_nPE'+nPE+endLabel+'.root')
    #fpathPred['pred_2017_corrih'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_corrTemplateIh_nPE'+nPE+endLabel+'.root')
    #fpathPred['pred_2017_corrmom'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_corrTemplateP_nPE'+nPE+endLabel+'.root')
    ### TODO fpathPred['pred_2017_fitihup'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_fitIhUp_nPE'+nPE+endLabel+'.root')
    ### TODO fpathPred['pred_2017_fitihdown'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_fitIhDown_nPE'+nPE+endLabel+'.root')
    ### TODO fpathPred['pred_2017_fitmomup'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_fitPUp_nPE'+nPE+endLabel+'.root')
    ### TODO fpathPred['pred_2017_fitmomdown'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_fitPDown_nPE'+nPE+endLabel+'.root')

    fpathPred['obs_2017'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta5_rebinIh4_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    #print(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta5_rebinIh4_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_nominal'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta5_rebinIh4_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_etaup'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_etadown'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta8_rebinIh4_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_ihup'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta5_rebinIh2_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_ihdown'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta5_rebinIh8_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_momup'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta5_rebinIh4_rebinP1_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_momdown'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta5_rebinIh4_rebinP4_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_corrih'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta5_rebinIh4_rebinP2_rebinMass1_corrTemplateIh_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_corrmom'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta5_rebinIh4_rebinP2_rebinMass1_corrTemplateP_nPE'+nPE+endLabel+'.root')
    #fpathPred['pred_2017_fitihup'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta5_rebinIh4_rebinP2_rebinMass1_fitIhUp_nPE'+nPE+endLabel+'.root')
    #fpathPred['pred_2017_fitihdown'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta5_rebinIh4_rebinP2_rebinMass1_fitIhDown_nPE'+nPE+endLabel+'.root')
    #fpathPred['pred_2017_fitmomup'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta5_rebinIh4_rebinP2_rebinMass1_fitPUp_nPE'+nPE+endLabel+'.root')
    #fpathPred['pred_2017_fitmomdown'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta5_rebinIh4_rebinP2_rebinMass1_fitPDown_nPE'+nPE+endLabel+'.root')

    year='2018'

    #fpathPred['pred_2018_nominal'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    #fpathPred['pred_2018_etaup'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta2_rebinIh4_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    #fpathPred['pred_2018_etadown'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta8_rebinIh4_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    #fpathPred['pred_2018_ihup'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh2_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    #fpathPred['pred_2018_ihdown'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh8_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    #fpathPred['pred_2018_momup'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP1_rebinMass1_nPE'+nPE+endLabel+'.root')
    #fpathPred['pred_2018_momdown'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP4_rebinMass1_nPE'+nPE+endLabel+'.root')
    #fpathPred['pred_2018_corrih'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_corrTemplateIh_nPE'+nPE+endLabel+'.root')
    #fpathPred['pred_2018_corrmom'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_corrTemplateP_nPE'+nPE+endLabel+'.root')
    ### TODO fpathPred['pred_2018_fitihup'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_fitIhUp_nPE'+nPE+endLabel+'.root')
    ### TODO fpathPred['pred_2018_fitihdown'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_fitIhDown_nPE'+nPE+endLabel+'.root')
    ### TODO fpathPred['pred_2018_fitmomdown'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_fitPDown_nPE'+nPE+endLabel+'.root')
    ### TODO fpathPred['pred_2018_fitmomup'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_fitPUp_nPE'+nPE+endLabel+'.root')
    ### TODO fpathPred['pred_2018_fitmomup2'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_fitPUp_nPE'+nPE+endLabel+'.root')
    
    fpathPred['obs_2018'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta5_rebinIh4_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_nominal'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta5_rebinIh4_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_etaup'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_etadown'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta8_rebinIh4_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_ihup'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta5_rebinIh2_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_ihdown'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta5_rebinIh8_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_momup'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta5_rebinIh4_rebinP1_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_momdown'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta5_rebinIh4_rebinP4_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_corrih'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta5_rebinIh4_rebinP2_rebinMass1_corrTemplateIh_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_corrmom'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta5_rebinIh4_rebinP2_rebinMass1_corrTemplateP_nPE'+nPE+endLabel+'.root')
    #fpathPred['pred_2018_fitihup'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta5_rebinIh4_rebinP2_rebinMass1_fitIhUp_nPE'+nPE+endLabel+'.root')
    #fpathPred['pred_2018_fitihdown'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta5_rebinIh4_rebinP2_rebinMass1_fitIhDown_nPE'+nPE+endLabel+'.root')
    #fpathPred['pred_2018_fitmomdown'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta5_rebinIh4_rebinP2_rebinMass1_fitPDown_nPE'+nPE+endLabel+'.root')
    #fpathPred['pred_2018_fitmomup'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta5_rebinIh4_rebinP2_rebinMass1_fitPUp_nPE'+nPE+endLabel+'.root')
    #fpathPred['pred_2018_fitmomup2'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta5_rebinIh4_rebinP2_rebinMass1_fitPUp_nPE'+nPE+endLabel+'.root')

    #2017 bias corr parameters
    if(searchRegion=='SR1'):
        p1=0.00108083
        p0=0.983476
    if(searchRegion=='SR2'):
        p1=0.00117178
        p0=0.981479
    if(searchRegion=='SR3'):
        #p1=0.000954989
        #p0=0.962746
        ### UNBLINDING V1
        p1=0.00131767
        p0=0.939438
        

    mass_plot['obs_2017'] = fpathPred['obs_2017'].Get("mass_obs_"+regionBckg)
    mass_plot['pred_2017_nominal'] = fpathPred['pred_2017_nominal'].Get("mass_predBC_"+regionBckg)
    mass_plot['pred_2017_stat'] = statErr(mass_plot['pred_2017_nominal'],"statErr")
    mass_plot['pred_2017_etaup'] = fpathPred['pred_2017_etaup'].Get("mass_predBC_"+regionBckg)
    mass_plot['pred_2017_etadown'] = fpathPred['pred_2017_etadown'].Get("mass_predBC_"+regionBckg)
    mass_plot['pred_2017_ihup'] = fpathPred['pred_2017_ihup'].Get("mass_predBC_"+regionBckg)
    mass_plot['pred_2017_ihdown'] = fpathPred['pred_2017_ihdown'].Get("mass_predBC_"+regionBckg)
    mass_plot['pred_2017_momup'] = fpathPred['pred_2017_momup'].Get("mass_predBC_"+regionBckg)
    mass_plot['pred_2017_momdown'] = fpathPred['pred_2017_momdown'].Get("mass_predBC_"+regionBckg)
    mass_plot['pred_2017_corrih'] = fpathPred['pred_2017_corrih'].Get("mass_predBC_"+regionBckg)
    mass_plot['pred_2017_corrmom'] = fpathPred['pred_2017_corrmom'].Get("mass_predBC_"+regionBckg)
    #mass_plot['pred_2017_fitihup'] = fpathPred['pred_2017_fitihup'].Get("mass_predBC_"+regionBckg)
    #mass_plot['pred_2017_fitihdown'] = fpathPred['pred_2017_fitihdown'].Get("mass_predBC_"+regionBckg)
    #mass_plot['pred_2017_fitmomup'] = fpathPred['pred_2017_fitmomup'].Get("mass_predBC_"+regionBckg)
    #mass_plot['pred_2017_fitmomdown'] = fpathPred['pred_2017_fitmomdown'].Get("mass_predBC_"+regionBckg)
    
    mass_plot['pred_2017_nominal'] = BiasCorrection(mass_plot['pred_2017_nominal'],p1,p0)
    mass_plot['pred_2017_stat'] = BiasCorrection(mass_plot['pred_2017_stat'],p1,p0)
    mass_plot['pred_2017_etaup'] = BiasCorrection(mass_plot['pred_2017_etaup'],p1,p0)
    mass_plot['pred_2017_etadown'] = BiasCorrection(mass_plot['pred_2017_etadown'],p1,p0)
    mass_plot['pred_2017_ihup'] = BiasCorrection(mass_plot['pred_2017_ihup'],p1,p0)
    mass_plot['pred_2017_ihdown'] = BiasCorrection(mass_plot['pred_2017_ihdown'],p1,p0)
    mass_plot['pred_2017_momup'] = BiasCorrection(mass_plot['pred_2017_momup'],p1,p0)
    mass_plot['pred_2017_momdown'] = BiasCorrection(mass_plot['pred_2017_momdown'],p1,p0)
    mass_plot['pred_2017_corrih'] = BiasCorrection(mass_plot['pred_2017_corrih'],p1,p0)
    mass_plot['pred_2017_corrmom'] = BiasCorrection(mass_plot['pred_2017_corrmom'],p1,p0)
    #mass_plot['pred_2017_fitihup'] = BiasCorrection(mass_plot['pred_2017_fitihup'],p1,p0)
    #mass_plot['pred_2017_fitihdown'] = BiasCorrection(mass_plot['pred_2017_fitihdown'],p1,p0)
    #mass_plot['pred_2017_fitmomup'] = BiasCorrection(mass_plot['pred_2017_fitmomup'],p1,p0)
    #mass_plot['pred_2017_fitmomdown'] = BiasCorrection(mass_plot['pred_2017_fitmomdown'],p1,p0)
    mass_plot['pred_2017_corrbias'] = BiasCorrection(mass_plot['pred_2017_nominal'],p1,p0)

    #2018 bias corr parameters
    if(searchRegion=='SR1'):
        p1=0.00121927
        p0=0.964028
    if(searchRegion=='SR2'):
        p1=0.00125541
        p0=0.96542
    if(searchRegion=='SR3'):
        #p1=0.00121957
        #p0=0.963581
        ### UNBLINDING V1
        p1=0.00169109
        p0=0.901472

    mass_plot['obs_2018'] = fpathPred['obs_2018'].Get("mass_obs_"+regionBckg)
    mass_plot['pred_2018_nominal'] = fpathPred['pred_2018_nominal'].Get("mass_predBC_"+regionBckg)
    mass_plot['pred_2018_etaup'] = fpathPred['pred_2018_etaup'].Get("mass_predBC_"+regionBckg)
    mass_plot['pred_2018_etadown'] = fpathPred['pred_2018_etadown'].Get("mass_predBC_"+regionBckg)
    mass_plot['pred_2018_ihup'] = fpathPred['pred_2018_ihup'].Get("mass_predBC_"+regionBckg)
    mass_plot['pred_2018_ihdown'] = fpathPred['pred_2018_ihdown'].Get("mass_predBC_"+regionBckg)
    mass_plot['pred_2018_momup'] = fpathPred['pred_2018_momup'].Get("mass_predBC_"+regionBckg)
    mass_plot['pred_2018_momdown'] = fpathPred['pred_2018_momdown'].Get("mass_predBC_"+regionBckg)
    mass_plot['pred_2018_corrih'] = fpathPred['pred_2018_corrih'].Get("mass_predBC_"+regionBckg)
    mass_plot['pred_2018_corrmom'] = fpathPred['pred_2018_corrmom'].Get("mass_predBC_"+regionBckg)
    #mass_plot['pred_2018_fitihup'] = fpathPred['pred_2018_fitihup'].Get("mass_predBC_"+regionBckg)
    #mass_plot['pred_2018_fitihdown'] = fpathPred['pred_2018_fitihdown'].Get("mass_predBC_"+regionBckg)
    #mass_plot['pred_2018_fitmomup'] = fpathPred['pred_2018_fitmomup'].Get("mass_predBC_"+regionBckg)
    #mass_plot['pred_2018_fitmomdown'] = fpathPred['pred_2018_fitmomdown'].Get("mass_predBC_"+regionBckg)
    #mass_plot['pred_2018_fitmomdown'] = fpathPred['pred_2018_fitmomdown'].Get("mass_predBC_"+regionBckg)
    mass_plot['pred_2018_stat'] = statErr(mass_plot['pred_2018_nominal'],"statErr")


    mass_plot['pred_2018_nominal'] = BiasCorrection(mass_plot['pred_2018_nominal'],p1,p0)
    mass_plot['pred_2018_stat'] = BiasCorrection(mass_plot['pred_2018_stat'],p1,p0)
    mass_plot['pred_2018_etaup'] = BiasCorrection(mass_plot['pred_2018_etaup'],p1,p0)
    mass_plot['pred_2018_etadown'] = BiasCorrection(mass_plot['pred_2018_etadown'],p1,p0)
    mass_plot['pred_2018_ihup'] = BiasCorrection(mass_plot['pred_2018_ihup'],p1,p0)
    mass_plot['pred_2018_ihdown'] = BiasCorrection(mass_plot['pred_2018_ihdown'],p1,p0)
    mass_plot['pred_2018_momup'] = BiasCorrection(mass_plot['pred_2018_momup'],p1,p0)
    mass_plot['pred_2018_momdown'] = BiasCorrection(mass_plot['pred_2018_momdown'],p1,p0)
    mass_plot['pred_2018_corrih'] = BiasCorrection(mass_plot['pred_2018_corrih'],p1,p0)
    mass_plot['pred_2018_corrmom'] = BiasCorrection(mass_plot['pred_2018_corrmom'],p1,p0)
    #mass_plot['pred_2018_fitihup'] = BiasCorrection(mass_plot['pred_2018_fitihup'],p1,p0)
    #mass_plot['pred_2018_fitihdown'] = BiasCorrection(mass_plot['pred_2018_fitihdown'],p1,p0)
    #mass_plot['pred_2018_fitmomup'] = BiasCorrection(mass_plot['pred_2018_fitmomup'],p1,p0)
    #mass_plot['pred_2018_fitmomdown'] = BiasCorrection(mass_plot['pred_2018_fitmomdown'],p1,p0)
    mass_plot['pred_2018_corrbias'] = BiasCorrection(mass_plot['pred_2018_nominal'],p1,p0)

    ofileBase.cd()
    #mass_plot['pred_2018_nominal'].SetName("nominal")
    #mass_plot['pred_2018_etaup'].SetName("EtaUp")
    #mass_plot['pred_2018_etadown'].SetName("EtaDown")
#
    #mass_plot['pred_2018_nominal'].Write()
    #mass_plot['pred_2018_etaup'].Write()
    #mass_plot['pred_2018_etadown'].Write()
    
    
    name = {
       'Gluino500_2018': '$\\tilde{g}$ (M=500 GeV)',
       'Gluino800_2018': '$\\tilde{g}$ (M=800 GeV)',
       'Gluino1000_2018': '$\\tilde{g}$ (M=1000 GeV)',
       'Gluino1400_2018': '$\\tilde{g}$ (M=1400 GeV)',
       'Gluino1600_2018': '$\\tilde{g}$ (M=1600 GeV)',
       'Gluino1800_2018': '$\\tilde{g}$ (M=1800 GeV)',
       'Gluino2000_2018': '$\\tilde{g}$ (M=2000 GeV)',
        'Gluino2200_2018': '$\\tilde{g}$ (M=2200 GeV)',
        'Gluino2400_2018': '$\\tilde{g}$ (M=2400 GeV)',
        'Gluino2600_2018': '$\\tilde{g}$ (M=2600 GeV)',
        'Stop800_2018': '$\\tilde{t}$ (M=800 GeV)',
        'Stop1000_2018': '$\\tilde{t}$ (M=1000 GeV)',
        'Stop1200_2018': '$\\tilde{t}$ (M=1200 GeV)',
       'Stop1400_2018': '$\\tilde{t}$ (M=1400 GeV)',
       'Stop1600_2018': '$\\tilde{t}$ (M=1600 GeV)',
       'Stop1800_2018': '$\\tilde{t}$ (M=1800 GeV)',
       'Stop2000_2018': '$\\tilde{t}$ (M=2000 GeV)',
        'Stop2200_2018': '$\\tilde{t}$ (M=2200 GeV)',
        'Stop2400_2018': '$\\tilde{t}$ (M=2400 GeV)',
        'Stop2600_2018': '$\\tilde{t}$ (M=2600 GeV)',
        'pairStau247_2018': '$\\tilde{\\tau}$ (M=247 GeV)',
        'pairStau308_2018': '$\\tilde{\\tau}$ (M=308 GeV)',
        'pairStau432_2018': '$\\tilde{\\tau}$ (M=432 GeV)',
        'pairStau557_2018': '$\\tilde{\\tau}$ (M=557 GeV)',
        'pairStau651_2018': '$\\tilde{\\tau}$ (M=651 GeV)',
        'pairStau745_2018': '$\\tilde{\\tau}$ (M=745 GeV)',
        'pairStau871_2018': '$\\tilde{\\tau}$ (M=871 GeV)',
        'pairStau1029_2018': '$\\tilde{\\tau}$ (M=1029 GeV)',
        'DYcharge1e_100_2018': '$\\tilde{\\tau}\'^{1e}$ (M=100 GeV)',
        'DYcharge1e_200_2018': '$\\tilde{\\tau}\'^{1e}$ (M=200 GeV)',
        'DYcharge1e_400_2018': '$\\tilde{\\tau}\'^{1e}$ (M=400 GeV)',
        'DYcharge1e_500_2018': '$\\tilde{\\tau}\'^{1e}$ (M=500 GeV)',
        'DYcharge1e_800_2018': '$\\tilde{\\tau}\'^{1e}$ (M=800 GeV)',
        'DYcharge1e_1000_2018': '$\\tilde{\\tau}\'^{1e}$ (M=1000 GeV)',
        'DYcharge1e_1400_2018': '$\\tilde{\\tau}\'^{1e}$ (M=1400 GeV)',
        'DYcharge1e_1800_2018': '$\\tilde{\\tau}\'^{1e}$ (M=1800 GeV)',
        'DYcharge1e_2200_2018': '$\\tilde{\\tau}\'^{1e}$ (M=2200 GeV)',
        'DYcharge1e_2600_2018': '$\\tilde{\\tau}\'^{1e}$ (M=2600 GeV)',
        'DYcharge2e_100_2018': '$\\tilde{\\tau}\'^{2e}$ (M=100 GeV)',
        'DYcharge2e_200_2018': '$\\tilde{\\tau}\'^{2e}$ (M=200 GeV)',
        'DYcharge2e_400_2018': '$\\tilde{\\tau}\'^{2e}$ (M=400 GeV)',
        'DYcharge2e_500_2018': '$\\tilde{\\tau}\'^{2e}$ (M=500 GeV)',
        'DYcharge2e_800_2018': '$\\tilde{\\tau}\'^{2e}$ (M=800 GeV)',
        'DYcharge2e_1000_2018': '$\\tilde{\\tau}\'^{2e}$ (M=1000 GeV)',
        'DYcharge2e_1400_2018': '$\\tilde{\\tau}\'^{2e}$ (M=1400 GeV)',
        'DYcharge2e_1800_2018': '$\\tilde{\\tau}\'^{2e}$ (M=1800 GeV)',
        'DYcharge2e_2200_2018': '$\\tilde{\\tau}\'^{2e}$ (M=2200 GeV)',
        'DYcharge2e_2600_2018': '$\\tilde{\\tau}\'^{2e}$ (M=2600 GeV)',

    }

    target={
       'Gluino500_2018': 500,
       'Gluino800_2018': 800,
       'Gluino1000_2018': 1000,
       'Gluino1400_2018': 1400,
       'Gluino1600_2018': 1600,
       'Gluino1800_2018': 1800,
       'Gluino2000_2018': 2000,
        'Gluino2200_2018': 2200,
        'Gluino2400_2018': 2400,
        'Gluino2600_2018': 2600,
       'Stop800_2018': 800,
       'Stop1000_2018': 1000,
       'Stop1200_2018': 1200,
       'Stop1400_2018': 1400,
       'Stop1600_2018': 1600,
       'Stop1800_2018': 1800,
       'Stop2000_2018': 2000,
        'Stop2200_2018': 2200,
        'Stop2400_2018': 2400,
        'Stop2600_2018': 2600,
        'pairStau308_2018': 308,
        'pairStau432_2018': 432,
        'pairStau557_2018': 557,
        'pairStau651_2018': 651,
        'pairStau745_2018': 745,
        'pairStau871_2018': 871,
        'pairStau1029_2018': 1029,
        'gmsbStau432_2018': 432,
        'gmsbStau557_2018': 557,
        'gmsbStau651_2018': 651,
        'gmsbStau745_2018': 745,
        'gmsbStau871_2018': 871,
        'gmsbStau1029_2018': 1029,
        'DYcharge1e_100_2018': 100,
        'DYcharge1e_200_2018': 200,
        'DYcharge1e_400_2018': 400,
        'DYcharge1e_500_2018': 500,
        'DYcharge1e_800_2018': 800,
        'DYcharge1e_1000_2018': 1000,
        'DYcharge1e_1400_2018': 1400,
        'DYcharge1e_1800_2018': 1800,
        'DYcharge1e_2200_2018': 2200,
        'DYcharge1e_2600_2018': 2600,
        'DYcharge2e_100_2018': 100,
        'DYcharge2e_200_2018': 200,
        'DYcharge2e_400_2018': 400,
        'DYcharge2e_500_2018': 500,
        'DYcharge2e_800_2018': 800,
        'DYcharge2e_1000_2018': 1000,
        'DYcharge2e_1400_2018': 1400,
        'DYcharge2e_1800_2018': 1800,
        'DYcharge2e_2200_2018': 2200,
        'DYcharge2e_2600_2018': 2600,

    }

    fpathPred['pred'] = '/opt/sbg/cms/safe1/cms/dapparu/HSCP/CMSSW_10_6_27/src/SUSYBSMAnalysis/BackgroundPrediction/predWithUpAndDown'+regionBckg+'_2017_2018.root'
    fPred = rt.TFile.Open(fpathPred['pred'], 'READ')
    mass_plot['predNominal'] = fPred.Get('pred_new')
    mass_plot['predUp'] = fPred.Get('predU')
    mass_plot['predDown'] = fPred.Get('predD')
    
    outDataCardsDir = "datacards_UnB_v1_"+searchRegion+"_Aug29_symmetricSyst/"
    os.system("mkdir -p {0}".format(outDataCardsDir))

    h2=rt.TH2F("h2",";Target Mass [GeV];Mass Window [GeV];[a.u.]",120,0,3000,400,0,4000)

    x = []
    yPU = []
    yFpix = []
    yGstrip = []
    yPt = []
    yTrigger = []
    yK = []
    yC = []
    yTotal = []

    for signal in fpath.keys():
        
        print signal
        
        ifile=rt.TFile(fpath[signal])
        nominalSignal=ifile.Get(massPlotsSignal['Nominal'])
        
        nominalSignal.SetName("signal_yield_"+signal)
        ofileBase.cd()
        nominalSignal.Write()
        
        mean=nominalSignal.GetMean()
        stddev=nominalSignal.GetStdDev()

        targetMass=target[signal]
        if("Stop" in signal):
            targetMass-=5

        fillH2(h2,targetMass,mean,stddev,signal)

        xmin=mean-stddev
        xmax=mean+2*stddev
        #xmax=100000

        if (xmin<300):
            xmin=300
        
        obs_2017 = integralHisto(mass_plot['obs_2017'],xmin,xmax)
        obs_2018 = integralHisto(mass_plot['obs_2018'],xmin,xmax)

        bkg_2017_nominal = integralHisto(mass_plot['pred_2017_nominal'],xmin,xmax)
        bkg_2017_stat = integralHisto(mass_plot['pred_2017_stat'],xmin,xmax)
        bkg_2017_etaBinning_up = integralHisto(mass_plot['pred_2017_etaup'],xmin,xmax)
        bkg_2017_etaBinning_down = integralHisto(mass_plot['pred_2017_etadown'],xmin,xmax)
        bkg_2017_ihBinning_up = integralHisto(mass_plot['pred_2017_ihup'],xmin,xmax)
        bkg_2017_ihBinning_down = integralHisto(mass_plot['pred_2017_ihdown'],xmin,xmax)
        bkg_2017_momBinning_up = integralHisto(mass_plot['pred_2017_momup'],xmin,xmax)
        bkg_2017_momBinning_down = integralHisto(mass_plot['pred_2017_momdown'],xmin,xmax)
        bkg_2017_corrTemplateIh_up = integralHisto(mass_plot['pred_2017_corrih'],xmin,xmax)
        bkg_2017_corrTemplateIh_down = integralHisto(mass_plot['pred_2017_corrih'],xmin,xmax)
        bkg_2017_corrTemplateMom_up = integralHisto(mass_plot['pred_2017_corrmom'],xmin,xmax)
        bkg_2017_corrTemplateMom_down = integralHisto(mass_plot['pred_2017_corrmom'],xmin,xmax)
        bkg_2017_fitIh_up = 0
        bkg_2017_fitIh_down = 0
        bkg_2017_fitMom_up = 0
        bkg_2017_fitMom_down = 0
        bkg_2017_correctionBias = integralHisto(mass_plot['pred_2017_corrbias'],xmin,xmax)

        bkg_2018_nominal = integralHisto(mass_plot['pred_2018_nominal'],xmin,xmax)
        bkg_2018_stat = integralHisto(mass_plot['pred_2018_stat'],xmin,xmax)
        bkg_2018_etaBinning_up = integralHisto(mass_plot['pred_2018_etaup'],xmin,xmax)
        bkg_2018_etaBinning_down = integralHisto(mass_plot['pred_2018_etadown'],xmin,xmax)
        bkg_2018_ihBinning_up = integralHisto(mass_plot['pred_2018_ihup'],xmin,xmax)
        bkg_2018_ihBinning_down = integralHisto(mass_plot['pred_2018_ihdown'],xmin,xmax)
        bkg_2018_momBinning_up = integralHisto(mass_plot['pred_2018_momup'],xmin,xmax)
        bkg_2018_momBinning_down = integralHisto(mass_plot['pred_2018_momdown'],xmin,xmax)
        bkg_2018_corrTemplateIh_up = integralHisto(mass_plot['pred_2018_corrih'],xmin,xmax)
        bkg_2018_corrTemplateIh_down = integralHisto(mass_plot['pred_2018_corrih'],xmin,xmax)
        bkg_2018_corrTemplateMom_up = integralHisto(mass_plot['pred_2018_corrmom'],xmin,xmax)
        bkg_2018_corrTemplateMom_down = integralHisto(mass_plot['pred_2018_corrmom'],xmin,xmax)
        bkg_2018_fitIh_up = 0
        bkg_2018_fitIh_down = 0
        bkg_2018_fitMom_up = 0
        bkg_2018_fitMom_down = 0
        bkg_2018_correctionBias = integralHisto(mass_plot['pred_2018_corrbias'],xmin,xmax)

        signal_yield=integralHisto(nominalSignal,xmin,xmax)
        signal_pu_up=integralHisto(ifile.Get(massPlotsSignal['PU_up']),xmin,xmax)
        signal_pu_down=integralHisto(ifile.Get(massPlotsSignal['PU_down']),xmin,xmax)
        signal_Fpix_up=integralHisto(ifile.Get(massPlotsSignal['Fpix_up']),xmin,xmax)
        signal_Fpix_down=integralHisto(ifile.Get(massPlotsSignal['Fpix_down']),xmin,xmax)
        signal_Gstrip_up=integralHisto(ifile.Get(massPlotsSignal['Gstrip_up']),xmin,xmax)
        signal_Gstrip_down=integralHisto(ifile.Get(massPlotsSignal['Gstrip_down']),xmin,xmax)
        signal_Pt_up=integralHisto(ifile.Get(massPlotsSignal['Pt_up']),xmin,xmax)
        signal_Pt_down=integralHisto(ifile.Get(massPlotsSignal['Pt_down']),xmin,xmax)
        signal_Trigger_up=integralHisto(ifile.Get(massPlotsSignal['Trigger_up']),xmin,xmax)
        signal_Trigger_down=integralHisto(ifile.Get(massPlotsSignal['Trigger_down']),xmin,xmax)
        signal_K_up=integralHisto(ifile.Get(massPlotsSignal['K_up']),xmin,xmax)
        signal_K_down=integralHisto(ifile.Get(massPlotsSignal['K_down']),xmin,xmax)
        signal_C_up=integralHisto(ifile.Get(massPlotsSignal['C_up']),xmin,xmax)
        signal_C_down=integralHisto(ifile.Get(massPlotsSignal['C_down']),xmin,xmax)

        yield_total=0
        yield_minus2sigma_minus1sigma=0
        yield_minus1sigma_mu=0
        yield_mu_plus1sigma=0
        yield_plus1sigma_plus2sigma=0
        yield_minus1sigma_plus2sigma=0

        print (int(mean-stddev), int(mean+2*stddev))
        print ("integral 2017: ", mass_plot['obs_2017'].Integral(mass_plot['obs_2017'].FindBin(mean-stddev),mass_plot['obs_2017'].FindBin(mean+2*stddev)))
        print ("integral 2018: ", mass_plot['obs_2018'].Integral(mass_plot['obs_2018'].FindBin(mean-stddev),mass_plot['obs_2018'].FindBin(mean+2*stddev)))

        '''
        for i in range(0,nominalSignal.GetNbinsX()+1):
            yield_total+=nominalSignal.GetBinContent(i)
            if (nominalSignal.GetBinLowEdge(i)>= mean-2*stddev and nominalSignal.GetBinLowEdge(i)<=mean-stddev):
                yield_minus2sigma_minus1sigma+=nominalSignal.GetBinContent(i)
            if (nominalSignal.GetBinLowEdge(i)>= mean-stddev and nominalSignal.GetBinLowEdge(i)<=mean):
                yield_minus1sigma_mu+=nominalSignal.GetBinContent(i)
            if (nominalSignal.GetBinLowEdge(i)>= mean and nominalSignal.GetBinLowEdge(i)<=mean+stddev):
                yield_mu_plus1sigma+=nominalSignal.GetBinContent(i)
            if (nominalSignal.GetBinLowEdge(i)>= mean+stddev and nominalSignal.GetBinLowEdge(i)<=mean+2*stddev):
                yield_plus1sigma_plus2sigma+=nominalSignal.GetBinContent(i)
            if (nominalSignal.GetBinLowEdge(i)>= mean-stddev and nominalSignal.GetBinLowEdge(i)<=mean+2*stddev):
                yield_minus1sigma_plus2sigma+=nominalSignal.GetBinContent(i)
                
            #if (nominalSignal.GetBinLowEdge(i)>= 1000 and nominalSignal.GetBinLowEdge(i)<=2500):
            if (nominalSignal.GetBinLowEdge(i)>= mean-stddev and nominalSignal.GetBinLowEdge(i)<=mean+2*stddev and nominalSignal.GetBinLowEdge(i)>=300):
            #if (nominalSignal.GetBinLowEdge(i)>= mean-2*stddev and nominalSignal.GetBinLowEdge(i)<=mean+2*stddev):
            #if(i>=nominalSignal.FindBin(mean-2*stddev) and i>=nominalSignal.FindBin(300)):
            
                signal_yield+=nominalSignal.GetBinContent(i)
                signal_pu_up+=ifile.Get(massPlotsSignal['PU_up']).GetBinContent(i)
                signal_pu_down+=ifile.Get(massPlotsSignal['PU_down']).GetBinContent(i)
                signal_Fpix_up+=ifile.Get(massPlotsSignal['Fpix_up']).GetBinContent(i)
                signal_Fpix_down+=ifile.Get(massPlotsSignal['Fpix_down']).GetBinContent(i)
                signal_Gstrip_up+=ifile.Get(massPlotsSignal['Gstrip_up']).GetBinContent(i)
                signal_Gstrip_down+=ifile.Get(massPlotsSignal['Gstrip_down']).GetBinContent(i)
                signal_Pt_up+=ifile.Get(massPlotsSignal['Pt_up']).GetBinContent(i)
                signal_Pt_down+=ifile.Get(massPlotsSignal['Pt_down']).GetBinContent(i)
                signal_Trigger_up+=ifile.Get(massPlotsSignal['Trigger_up']).GetBinContent(i)
                signal_Trigger_down+=ifile.Get(massPlotsSignal['Trigger_down']).GetBinContent(i)
                signal_K_up+=ifile.Get(massPlotsSignal['K_up']).GetBinContent(i)
                signal_K_down+=ifile.Get(massPlotsSignal['K_down']).GetBinContent(i)
                signal_C_up+=ifile.Get(massPlotsSignal['C_up']).GetBinContent(i)
                signal_C_down+=ifile.Get(massPlotsSignal['C_down']).GetBinContent(i)
        for i in range(0,mass_plot['predNominal'].GetNbinsX()+1):
            #print mass_plot['predNominal'].GetBinLowEdge(i), mean-stddev, mean+2*stddev
            if (mass_plot['predNominal'].GetBinLowEdge(i)>= mean-stddev and mass_plot['predNominal'].GetBinLowEdge(i)<=mean+2*stddev):
            #if(i>=mass_plot['predNominal'].FindBin(mean) and i<=mass_plot['predNominal'].FindBin(mean+2*stddev) and i>=mass_plot['predNominal'].FindBin(300)):
                bkg_nominal+=mass_plot['predNominal'].GetBinContent(i)
                bkg_up+=mass_plot['predUp'].GetBinContent(i)
                bkg_down+=mass_plot['predDown'].GetBinContent(i)
                #print bkg_nominal
        
        for i in range(0,mass_plot['pred_2017_nominal'].GetNbinsX()+1):
            #if(i>=mass_plot['pred_2017_nominal'].FindBin(1000) and i<=mass_plot['pred_2017_nominal'].FindBin(2500)):
            #print 'bin' , i, mean-stddev, mass_plot['pred_2017_nominal'].FindBin(mean-stddev), mass_plot['pred_2017_nominal'].FindBin(300), mass_plot['pred_2017_nominal'].GetBinLowEdge(i)
            #if(i>=mass_plot['obs_2017'].FindBin(300)):
            #    print mass_plot['obs_2017'].GetBinContent(i)
            if(i>=mass_plot['pred_2017_nominal'].FindBin(mean-stddev) and i<=mass_plot['pred_2017_nominal'].FindBin(mean+2*stddev) and i>=mass_plot['pred_2017_nominal'].FindBin(300)):
            #if(i>=mass_plot['pred_2017_nominal'].FindBin(mean-2*stddev) and i>=mass_plot['pred_2017_nominal'].FindBin(300)):
                #obs_2017+=mass_plot['obs_2017'].GetBinContent(i)
                print mass_plot['obs_2017'].GetBinContent(i), i, mass_plot['obs_2017'].GetBinLowEdge(i)
                bkg_2017_nominal+=mass_plot['pred_2017_nominal'].GetBinContent(i)
                bkg_2017_stat+=mass_plot['pred_2017_stat'].GetBinContent(i)
                bkg_2017_etaBinning_up+=mass_plot['pred_2017_etaup'].GetBinContent(i)
                bkg_2017_etaBinning_down+=mass_plot['pred_2017_etadown'].GetBinContent(i)
                bkg_2017_ihBinning_up+=mass_plot['pred_2017_ihup'].GetBinContent(i)
                bkg_2017_ihBinning_down+=mass_plot['pred_2017_ihdown'].GetBinContent(i)
                bkg_2017_momBinning_up+=mass_plot['pred_2017_momup'].GetBinContent(i)
                bkg_2017_momBinning_down+=mass_plot['pred_2017_momdown'].GetBinContent(i)
                bkg_2017_corrTemplateIh_up+=mass_plot['pred_2017_corrih'].GetBinContent(i)
                bkg_2017_corrTemplateIh_down+=mass_plot['pred_2017_corrih'].GetBinContent(i)
                bkg_2017_corrTemplateMom_up+=mass_plot['pred_2017_corrmom'].GetBinContent(i)
                bkg_2017_corrTemplateMom_down+=mass_plot['pred_2017_corrmom'].GetBinContent(i)
                #bkg_2017_fitIh_up+=mass_plot['pred_2017_fitihup'].GetBinContent(i)
                #bkg_2017_fitIh_down+=mass_plot['pred_2017_fitihdown'].GetBinContent(i)
                #bkg_2017_fitMom_up+=mass_plot['pred_2017_fitmomup'].GetBinContent(i)
                #bkg_2017_fitMom_down+=mass_plot['pred_2017_fitmomdown'].GetBinContent(i)
                bkg_2017_correctionBias+=mass_plot['pred_2017_corrbias'].GetBinContent(i)
                #print i,mass_plot['pred_2017_corrbias'].GetBinLowEdge(i),mass_plot['pred_2017_corrbias'].GetBinContent(i), mass_plot['pred_2017_nominal'].GetBinContent(i)*(0.000954989*mass_plot['pred_2017_nominal'].GetBinCenter(i)+0.962746),mass_plot['pred_2017_nominal'].GetBinContent(i),mass_plot['pred_2017_corrbias'].GetBinContent(i)/mass_plot['pred_2017_nominal'].GetBinContent(i)
                #print 'eta binning', bkg_2017_etaBinning_up, bkg_2017_etaBinning_down, bkg_2017_nominal, bkg_2017_etaBinning_up/bkg_2017_nominal, bkg_2017_etaBinning_down/bkg_2017_nominal
        for i in range(0,mass_plot['pred_2018_nominal'].GetNbinsX()+1):
            #if(i>=mass_plot['pred_2018_nominal'].FindBin(1000) and i<=mass_plot['pred_2018_nominal'].FindBin(2500)):
            if(i>=mass_plot['pred_2018_nominal'].FindBin(mean-stddev) and i<=mass_plot['pred_2018_nominal'].FindBin(mean+2*stddev) and i>=mass_plot['pred_2018_nominal'].FindBin(300)):
            #if(i>=mass_plot['pred_2018_nominal'].FindBin(mean-2*stddev) and i>=mass_plot['pred_2018_nominal'].FindBin(300)):
                
                #obs_2018+=mass_plot['obs_2018'].GetBinContent(i)
                bkg_2018_nominal+=mass_plot['pred_2018_nominal'].GetBinContent(i)
                bkg_2018_stat+=mass_plot['pred_2018_stat'].GetBinContent(i)
                bkg_2018_etaBinning_up+=mass_plot['pred_2018_etaup'].GetBinContent(i)
                bkg_2018_etaBinning_down+=mass_plot['pred_2018_etadown'].GetBinContent(i)
                bkg_2018_ihBinning_up+=mass_plot['pred_2018_ihup'].GetBinContent(i)
                bkg_2018_ihBinning_down+=mass_plot['pred_2018_ihdown'].GetBinContent(i)
                bkg_2018_momBinning_up+=mass_plot['pred_2018_momup'].GetBinContent(i)
                bkg_2018_momBinning_down+=mass_plot['pred_2018_momdown'].GetBinContent(i)
                bkg_2018_corrTemplateIh_up+=mass_plot['pred_2018_corrih'].GetBinContent(i)
                bkg_2018_corrTemplateIh_down+=mass_plot['pred_2018_corrih'].GetBinContent(i)
                bkg_2018_corrTemplateMom_up+=mass_plot['pred_2018_corrmom'].GetBinContent(i)
                bkg_2018_corrTemplateMom_down+=mass_plot['pred_2018_corrmom'].GetBinContent(i)
                #bkg_2018_fitIh_up+=mass_plot['pred_2018_fitihup'].GetBinContent(i)
                #bkg_2018_fitIh_down+=mass_plot['pred_2018_fitihdown'].GetBinContent(i)
                #bkg_2018_fitMom_up+=mass_plot['pred_2018_fitmomup'].GetBinContent(i)
                #bkg_2018_fitMom_down+=mass_plot['pred_2018_fitmomdown'].GetBinContent(i)
                bkg_2018_correctionBias+=mass_plot['pred_2018_corrbias'].GetBinContent(i)
        '''
        #print '2017', mass_plot['pred_2017_fitmomdown'].Integral(), mass_plot['pred_2017_fitmomup'].Integral(),mass_plot['pred_2017_nominal'].Integral()
        #print '2018', mass_plot['pred_2018_fitmomdown'].Integral(), mass_plot['pred_2018_fitmomup'].Integral(),mass_plot['pred_2018_nominal'].Integral()
        
        #print yield_total/yield_total, yield_minus2sigma_minus1sigma/yield_total, yield_minus1sigma_mu/yield_total, yield_mu_plus1sigma/yield_total, yield_plus1sigma_plus2sigma/yield_total, yield_minus1sigma_plus2sigma/yield_total
        
        #obs_2017 = mass_plot['obs_2017'].Integral(mass_plot['obs_2017'].FindBin(mean-stddev),mass_plot['obs_2017'].FindBin(mean+2*stddev))
        #obs_2018 = mass_plot['obs_2018'].Integral(mass_plot['obs_2018'].FindBin(mean-stddev),mass_plot['obs_2018'].FindBin(mean+2*stddev))

        #if(mean-stddev<300):
        #    obs_2017 = mass_plot['obs_2017'].Integral(mass_plot['obs_2017'].FindBin(300),mass_plot['obs_2017'].FindBin(mean+2*stddev))
        #    obs_2018 = mass_plot['obs_2018'].Integral(mass_plot['obs_2018'].FindBin(300),mass_plot['obs_2018'].FindBin(mean+2*stddev))

        if(signal_yield!=0):
            signal_pu_up/=signal_yield
            signal_pu_down/=signal_yield
            signal_Fpix_up/=signal_yield
            signal_Fpix_down/=signal_yield
            signal_Gstrip_up/=signal_yield
            signal_Gstrip_down/=signal_yield
            signal_Pt_up/=signal_yield
            signal_Pt_down/=signal_yield
            signal_Trigger_up/=signal_yield
            signal_Trigger_down/=signal_yield
            signal_K_up/=signal_yield
            signal_K_down/=signal_yield
            signal_C_up/=signal_yield
            signal_C_down/=signal_yield

        if(systSignal in signal):

            x.append(targetMass)

            yPU.append(pushSyst(signal_pu_up,signal_pu_down))
            yFpix.append(pushSyst(signal_Fpix_up,signal_Fpix_down))
            yGstrip.append(pushSyst(signal_Gstrip_up,signal_Gstrip_down))
            yPt.append(pushSyst(signal_Pt_up,signal_Pt_down))
            yTrigger.append(pushSyst(signal_Trigger_up,signal_Trigger_down))
            yK.append(pushSyst(signal_K_up,signal_K_down))
            yC.append(pushSyst(signal_C_up,signal_C_down))

            yTotal.append(np.sqrt(pow(pushSyst(signal_pu_up,signal_pu_down),2)+pow(pushSyst(signal_Fpix_up,signal_Fpix_down),2)+pow(pushSyst(signal_Gstrip_up,signal_Gstrip_down),2)+pow(pushSyst(signal_Pt_up,signal_Pt_down),2)+pow(pushSyst(signal_Trigger_up,signal_Trigger_down),2)+pow(pushSyst(signal_K_up,signal_K_down),2)+pow(pushSyst(signal_C_up,signal_C_down),2)))


        if(bkg_2017_nominal!=0):
            bkg_2017_stat/=bkg_2017_nominal
            bkg_2017_etaBinning_up/=bkg_2017_nominal
            bkg_2017_etaBinning_down/=bkg_2017_nominal
            bkg_2017_ihBinning_up/=bkg_2017_nominal
            bkg_2017_ihBinning_down/=bkg_2017_nominal
            bkg_2017_momBinning_up/=bkg_2017_nominal
            bkg_2017_momBinning_down/=bkg_2017_nominal
            bkg_2017_corrTemplateIh_up/=bkg_2017_nominal
            bkg_2017_corrTemplateIh_down/=bkg_2017_nominal
            bkg_2017_corrTemplateMom_up/=bkg_2017_nominal
            bkg_2017_corrTemplateMom_down/=bkg_2017_nominal
            #bkg_2017_fitIh_up/=bkg_2017_nominal
            #bkg_2017_fitIh_down/=bkg_2017_nominal
            #bkg_2017_fitMom_up/=bkg_2017_nominal
            #bkg_2017_fitMom_down/=bkg_2017_nominal
            bkg_2017_fitIh_up=1
            bkg_2017_fitIh_down=1
            bkg_2017_fitMom_up=1
            bkg_2017_fitMom_down=1
            bkg_2017_correctionBias/=bkg_2017_nominal

        
        if(bkg_2018_nominal!=0):
            bkg_2018_stat/=bkg_2018_nominal
            bkg_2018_etaBinning_up/=bkg_2018_nominal
            bkg_2018_etaBinning_down/=bkg_2018_nominal
            bkg_2018_ihBinning_up/=bkg_2018_nominal
            bkg_2018_ihBinning_down/=bkg_2018_nominal
            bkg_2018_momBinning_up/=bkg_2018_nominal
            bkg_2018_momBinning_down/=bkg_2018_nominal
            bkg_2018_corrTemplateIh_up/=bkg_2018_nominal
            bkg_2018_corrTemplateIh_down/=bkg_2018_nominal
            bkg_2018_corrTemplateMom_up/=bkg_2018_nominal
            bkg_2018_corrTemplateMom_down/=bkg_2018_nominal
            #bkg_2018_fitIh_up/=bkg_2018_nominal
            #bkg_2018_fitIh_down/=bkg_2018_nominal
            #bkg_2018_fitMom_up/=bkg_2018_nominal
            #bkg_2018_fitMom_down/=bkg_2018_nominal
            bkg_2018_fitIh_up=1
            bkg_2018_fitIh_down=1
            bkg_2018_fitMom_up=1
            bkg_2018_fitMom_down=1
            bkg_2018_correctionBias/=bkg_2018_nominal

        #print bkg_2018_fitMom_up,bkg_2018_fitMom_down

        sig_2017_unc = {
            'sig2017_Fpix': [signal_Fpix_down, signal_Fpix_up],
            'sig2017_Gstrip': [signal_Gstrip_down, signal_Gstrip_up],
            'sig2017_Pt': [signal_Pt_down, signal_Pt_up],
            'sig2017_K': [signal_K_down, signal_K_up],
            'sig2017_C': [signal_C_down, signal_C_up],
        }

        sig_2018_unc = {
            'sig2018_Fpix': [signal_Fpix_down, signal_Fpix_up],
            'sig2018_Gstrip': [signal_Gstrip_down, signal_Gstrip_up],
            'sig2018_Pt': [signal_Pt_down, signal_Pt_up],
            'sig2018_K': [signal_K_down, signal_K_up],
            'sig2018_C': [signal_C_down, signal_C_up],
        }

        #sig_2017_unc = {
        #    'sig2017_Fpix': [1],
        #    'sig2017_Gstrip': [1],
        #    'sig2017_Pt': [1],
        #    'sig2017_K': [1],
        #    'sig2017_C': [1],
        #}
#
        #sig_2018_unc = {
        #    'sig2018_Fpix': [1],
        #    'sig2018_Gstrip': [1],
        #    'sig2018_Pt': [1],
        #    'sig2018_K': [1],
        #    'sig2018_C': [1],
        #}

        sig_unc_correlated = {
            'sig_pu' : [signal_pu_down, signal_pu_up, signal_pu_down, signal_pu_up],
            'sig_Trigger' : [signal_Trigger_down, signal_Trigger_up, signal_Trigger_down, signal_Trigger_up],
        }

        #sig_unc_correlated = {
        #    'sig_pu' : [1,1,1,1],
        #    'sig_Trigger' : [1,1,1,1],
        #}

        tot_unc_sig_2017 = { 
            'sig2017_Fpix': [signal_Fpix_down, signal_Fpix_up],
            'sig2017_Gstrip': [signal_Gstrip_down, signal_Gstrip_up],
            'sig2017_Pt': [signal_Pt_down, signal_Pt_up],
            'sig2017_K': [signal_K_down, signal_K_up],
            'sig2017_C': [signal_C_down, signal_C_up],
            'sig_pu' : [signal_pu_down, signal_pu_up],
            'sig_Trigger' : [signal_Trigger_down, signal_Trigger_up],
        }

        tot_unc_sig_2018 = { 
            'sig2018_Fpix': [signal_Fpix_down, signal_Fpix_up],
            'sig2018_Gstrip': [signal_Gstrip_down, signal_Gstrip_up],
            'sig2018_Pt': [signal_Pt_down, signal_Pt_up],
            'sig2018_K': [signal_K_down, signal_K_up],
            'sig2018_C': [signal_C_down, signal_C_up],
            'sig_pu' : [signal_pu_down, signal_pu_up],
            'sig_Trigger' : [signal_Trigger_down, signal_Trigger_up],
        }

        bkg_2017_unc = {
            'bkg2017_stat': [bkg_2017_stat],
            'bkg2017_fitIh': [bkg_2017_fitIh_down, bkg_2017_fitIh_up],
            'bkg2017_fitMom': [bkg_2017_fitMom_down, bkg_2017_fitMom_up],
        }

        bkg_2018_unc = {
            'bkg2018_stat': [bkg_2018_stat],
            'bkg2018_fitIh': [bkg_2018_fitIh_down, bkg_2018_fitIh_up],
            'bkg2018_fitMom': [bkg_2018_fitMom_down, bkg_2018_fitMom_up],
        }

        #bkg_2017_unc = {
        #    'bkg2017_stat': [1.4],
        #    'bkg2017_fitIh': [1],
        #    'bkg2017_fitMom': [1],
        #}
#
        #bkg_2018_unc = {
        #    'bkg2018_stat': [1.4],
        #    'bkg2018_fitIh': [1],
        #    'bkg2018_fitMom': [1],
        #}

        bkg_unc_correlated = {
            'bkg_etaBinning' : [bkg_2017_etaBinning_down, bkg_2017_etaBinning_up, bkg_2018_etaBinning_down, bkg_2018_etaBinning_up],
            'bkg_ihBinning' : [bkg_2017_ihBinning_down, bkg_2017_ihBinning_up, bkg_2018_ihBinning_down, bkg_2018_ihBinning_up],
            'bkg_momentumBinning' : [bkg_2017_momBinning_down, bkg_2017_momBinning_up, bkg_2018_momBinning_down, bkg_2018_momBinning_up],
            'bkg_correctionTemplateIh' : [bkg_2017_corrTemplateIh_up, bkg_2018_corrTemplateIh_up],
            'bkg_correctionTemplateMomentum' : [bkg_2017_corrTemplateMom_up, bkg_2018_corrTemplateMom_up],
            'bkg_correctionBias' : [bkg_2017_correctionBias, bkg_2018_correctionBias],
        }
        
        #bkg_unc_correlated = {
        #    'bkg_etaBinning' : [1,1,1,1],
        #    'bkg_ihBinning' : [1,1,1,1],
        #    'bkg_momentumBinning' : [1,1,1,1],
        #    'bkg_correctionTemplateIh' : [1,1],
        #    'bkg_correctionTemplateMomentum' : [1,1],
        #    'bkg_correctionBias' : [1,1],
        #}

        tot_unc_bkg_2017 = {
            'bkg2017_stat': [bkg_2017_stat],
            'bkg2017_fitIh': [bkg_2017_fitIh_down, bkg_2017_fitIh_up],
            'bkg2017_fitMom': [bkg_2017_fitMom_down, bkg_2017_fitMom_up],
            'bkg_etaBinning' : [bkg_2017_etaBinning_down, bkg_2017_etaBinning_up],
            'bkg_ihBinning' : [bkg_2017_ihBinning_down, bkg_2017_ihBinning_up],
            'bkg_momentumBinning' : [bkg_2017_momBinning_down, bkg_2017_momBinning_up],
            'bkg_correctionTemplateIh' : [bkg_2017_corrTemplateIh_up],
            'bkg_correctionTemplateMomentum' : [bkg_2017_corrTemplateMom_up],
            'bkg_correctionBias' : [bkg_2017_correctionBias],
        }

        tot_unc_bkg_2018 = {
            'bkg2018_stat': [bkg_2018_stat],
            'bkg2018_fitIh': [bkg_2018_fitIh_down, bkg_2018_fitIh_up],
            'bkg2018_fitMom': [bkg_2018_fitMom_down, bkg_2018_fitMom_up],
            'bkg_etaBinning' : [bkg_2018_etaBinning_down, bkg_2018_etaBinning_up],
            'bkg_ihBinning' : [bkg_2018_ihBinning_down, bkg_2018_ihBinning_up],
            'bkg_momentumBinning' : [bkg_2018_momBinning_down, bkg_2018_momBinning_up],
            'bkg_correctionTemplateIh' : [bkg_2018_corrTemplateIh_up],
            'bkg_correctionTemplateMomentum' : [bkg_2018_corrTemplateMom_up],
            'bkg_correctionBias' : [bkg_2018_correctionBias],
        }

        bkg=bkg_2018_nominal
        bckg_up=totalUncertainy(tot_unc_bkg_2018,1)
        bckg_down=totalUncertainy(tot_unc_bkg_2018,0)
        signal_up=totalUncertainy(tot_unc_sig_2018,1)
        signal_down=totalUncertainy(tot_unc_sig_2018,0)
        #print bkg, bkg*bckg_down, bkg*bckg_up, signal_yield, signal_up, signal_down

        #signal_yield*= integralF


        signal_yield_norm=signal_yield*41.5/101
        signal_yield_norm=signal_yield*59.7/101
        #bckg_up*=bkg
        #bckg_down*=bkg
        #signal_up*=signal_yield_norm
        #signal_down*=signal_yield_norm



        make_datacard_hscp_combining2017and2018(outDataCardsDir,  signal, signal_yield*41.5/101, signal_yield*59.7/101, bkg_2017_nominal, bkg_2018_nominal, obs_2017, obs_2018, sig_2017_unc, sig_2018_unc, sig_unc_correlated, bkg_2017_unc, bkg_2018_unc, bkg_unc_correlated)
        #make_datacard_hscp_combining2017and2018(outDataCardsDir,  signal, signal_yield*41.5/101, signal_yield*59.7/101, bkg_2017_nominal, bkg_2018_nominal, bkg_2017_nominal, bkg_2018_nominal, sig_2017_unc, sig_2018_unc, sig_unc_correlated, bkg_2017_unc, bkg_2018_unc, bkg_unc_correlated)
        
        bkg=bkg_2017_nominal
        bckg_up=totalUncertainy(tot_unc_bkg_2017,1)
        bckg_down=totalUncertainy(tot_unc_bkg_2017,0)
        signal_up=totalUncertainy(tot_unc_sig_2017,1)
        signal_down=totalUncertainy(tot_unc_sig_2017,0)
        signal_yield_norm=signal_yield*41.5/101
        
        makeYieldFile(text_file_tex,name[signal],bkg,bckg_up,bckg_down,signal_yield_norm,signal_up,signal_down)
        
        bkg=bkg_2018_nominal
        bckg_up=totalUncertainy(tot_unc_bkg_2018,1)
        bckg_down=totalUncertainy(tot_unc_bkg_2018,0)
        signal_up=totalUncertainy(tot_unc_sig_2018,1)
        signal_down=totalUncertainy(tot_unc_sig_2018,0)
        signal_yield_norm=signal_yield*59.7/101
        
        makeYieldFile(text_file_tex_2018,name[signal],bkg,bckg_up,bckg_down,signal_yield_norm,signal_up,signal_down)

    grPU = create_TGraph(x, yPU, axis_title=['target mass [GeV]', 'Systematics uncertainties [%]'])
    grFpix = create_TGraph(x, yFpix, axis_title=['target mass [GeV]', 'Systematics uncertainties [%]'])
    grGstrip = create_TGraph(x, yGstrip, axis_title=['target mass [GeV]', 'Systematics uncertainties [%]'])
    grPt = create_TGraph(x, yPt, axis_title=['target mass [GeV]', 'Systematics uncertainties [%]'])
    grTrigger = create_TGraph(x, yTrigger, axis_title=['target mass [GeV]', 'Systematics uncertainties [%]'])
    grK = create_TGraph(x, yK, axis_title=['target mass [GeV]', 'Systematics uncertainties [%]'])
    grC = create_TGraph(x, yC, axis_title=['target mass [GeV]', 'Systematics uncertainties [%]'])
    grTotal = create_TGraph(x, yTotal, axis_title=['target mass [GeV]', 'Systematics uncertainties [%]'])


    text_file_tex.write('\n \\end{tabular}')
    text_file_tex.write('\n \\end{center}')
    text_file_tex.write('\n \\end{document}')

    text_file_tex.close()
            

    text_file_tex_2018.write('\n \\end{tabular}')
    text_file_tex_2018.write('\n \\end{center}')
    text_file_tex_2018.write('\n \\end{document}')

    text_file_tex_2018.close()
            

    tdrstyle.setTDRStyle()
    CMS_lumi.cmsText     = "CMS"
    iPos = 0
    CMS_lumi.extraText = "Internal"
    CMS_lumi.writeExtraText=True

    if( iPos==0 ): CMS_lumi.relPosX = 0.12
    # CMS_lumi.CMS_lumi(c, 4, 0)
    CMS_lumi.lumi_13TeV  = "101 fb^{-1}"
    
        
    rt.gStyle.SetOptStat(0)
    c1=rt.TCanvas()
    h2.Draw("col")
    h2.GetXaxis().SetTitleSize(0.04)
    h2.GetYaxis().SetTitleSize(0.04)
    h2.GetXaxis().SetTitleOffset(1.5)
    h2.GetYaxis().SetTitleOffset(2)
    CMS_lumi.CMS_lumi(c1, 4, iPos)
    c1.SaveAs("mass_window_Aug29/h2_massWindow"+regionSignal+".root")
    c1.SaveAs("mass_window_Aug29/h2_massWindow"+regionSignal+".pdf")

    setColorAndMarkerGr(grK,30,21)
    setColorAndMarkerGr(grC,38,22)
    setColorAndMarkerGr(grPU,46,23)
    setColorAndMarkerGr(grFpix,43,43)
    setColorAndMarkerGr(grGstrip,45,45)
    setColorAndMarkerGr(grPt,39,29)
    setColorAndMarkerGr(grTrigger,40,39)
    setColorAndMarkerGr(grTotal,28,34)

    leg2=rt.TLegend(0.5,0.3,0.7,0.5)
    leg2.AddEntry(grTotal,"Total","PE1")
    leg2.AddEntry(grK,"K","PE1")
    leg2.AddEntry(grC,"C","PE1")
    leg2.AddEntry(grPU,"PU","PE1")
    leg2.AddEntry(grFpix,"F^{pixel}","PE1")
    leg2.AddEntry(grGstrip,"G^{strip}","PE1")
    leg2.AddEntry(grPt,"p_{T}","PE1")
    leg2.AddEntry(grTrigger,"Trigger","PE1")

    c2=rt.TCanvas("c","c",800,800)
    c2.SetGridx()
    c2.SetGridy()
    grPU.Draw("AP")
    grPU.SetMinimum(0)
    grPU.SetMaximum(20)
    grPU.GetXaxis().SetLabelSize(0.03)
    grPU.GetXaxis().SetLabelOffset(0.015)
    grPU.GetXaxis().SetTitleSize(0.04)
    grPU.GetYaxis().SetTitleSize(0.04)
    grPU.GetXaxis().SetTitleOffset(1.5)
    grPU.GetYaxis().SetTitleOffset(2)
    grFpix.Draw("P")
    grGstrip.Draw("P")
    grPt.Draw("P")
    grTrigger.Draw("P")
    grK.Draw("P")
    grC.Draw("P")
    grTotal.Draw("P")
    leg2.Draw("same")
    CMS_lumi.CMS_lumi(c2, 4, iPos)
    c2.SaveAs("systTargetMass_Aug22/systTargetMass_"+systSignal+"_"+regionSignal+".root")
    c2.SaveAs("systTargetMass_Aug22/systTargetMass_"+systSignal+"_"+regionSignal+".pdf")
