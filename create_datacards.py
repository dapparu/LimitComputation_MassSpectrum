import ROOT as rt
import csv
import re
import sys
import collections
import os
sys.argv.append(' -b- ')

from collections import OrderedDict
import numpy as np
import array

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
            statErr.SetBinContent(i,h1.GetBinError(i)+h1.GetBinContent(i))
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

    #### uncertainties ####
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

def make_datacard_hscp_combining2017and2018_massShape(outDataCardsDir,  modelName, signal2017, signal2018, bkg2017, bkg2018, observation2017, observation2018, sig_2017_unc, sig_2018_unc, sig_unc_correlated, bkg_2017_unc, bkg_2018_unc, bkg_unc_correlated):

    text_file = open(outDataCardsDir+modelName+".txt", "w")

    text_file.write('imax {0} \n'.format(2))
    text_file.write('jmax {0} \n'.format(1))
    text_file.write('kmax * \n')
    text_file.write('shapes * * $PROCESS_$CHANNEL $PROCESS_$CHANNEL_$SYSTEMATIC \n')
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

    #### uncertainties ####
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


    text_file.close()

def totalUncertainy(unc,i=1):
    total=0
    for k,v in unc.items():
        if len(v)==2 and (i==0 or i==1):
            total+=pow(1-v[i],2)
        else:total+=pow(1-v[0],2)
    return np.sqrt(total)

def makeYieldFile(text_file_tex,modelName,bkg,bckg_up,bckg_down,obs,signal,signal_up,signal_down):
    typeOfDisplay='.2E'
    text_file_tex.write('\n '+modelName+' & $'+str(format(bkg,typeOfDisplay))+'^{+'+str(format(bkg*bckg_up,typeOfDisplay))+'}_{-'+str(format(bkg*bckg_down,typeOfDisplay))+'}$ & '+str(format(obs,typeOfDisplay))+' & $'+str(format(signal,typeOfDisplay))+'^{+'+str(format(signal*signal_up,typeOfDisplay))+'}_{-'+str(format(signal*signal_down,typeOfDisplay))+'}$ \\\\')
    text_file_tex.write('\n \hline')  
    
def fillH2(h2,targetMass,mean,stddev,s):
    xmin=mean-stddev
    if (xmin<300): 
        xmin=300
    xmax=mean+2*stddev
    for i in range(0,h2.GetNbinsX()+1):
        if (i!=h2.GetXaxis().FindBin(targetMass)): 
                continue
        else:
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

def saveSignalShape(outfile, h_signal, label):
    label = label.find("_")[0]
    print label
    h_2017 = h_signal.Clone("signal_Ch2017_"+label)
    h_2018 = h_signal.Clone("signal_Ch2018_"+label)
    if(label=="Nominal"): 
        h_2017 = h_signal.Clone("signal_Ch2017")
        h_2018 = h_signal.Clone("signal_Ch2018")
    h_2017.Scale(41.5/101.)
    h_2018.Scale(59.7/101.)
    outfile.cd()
    h_2017.Write()
    h_2018.Write()

def saveSigShape(outfile, h_signal, label):
    rebinning=array.array('d',[0.,20.,40.,60.,80.,100.,120.,140.,160.,180.,200.,220.,240.,260.,280.,300.,320.,340.,360.,380.,405.,435.,475.,525.,585.,660.,755.,875.,1025.,1210.,1440.,1730.,2000.,2500.,3200.,4000.])
    h_2017 = h_signal.Clone("signal_Ch2017_"+label)
    h_2018 = h_signal.Clone("signal_Ch2018_"+label)
    if(label=="Nominal"): 
        h_2017 = h_signal.Clone("signal_Ch2017")
        h_2018 = h_signal.Clone("signal_Ch2018")
    h_2017.Scale(41.5/101.)
    h_2018.Scale(59.7/101.)
    h_2017=h_2017.Rebin(len(rebinning)-1,"",rebinning)
    h_2018=h_2018.Rebin(len(rebinning)-1,"",rebinning)
    outfile.cd()
    h_2017.Write()
    h_2018.Write()

def savePredShape(outfile, h_pred, label):
    rebinning=array.array('d',[0.,20.,40.,60.,80.,100.,120.,140.,160.,180.,200.,220.,240.,260.,280.,300.,320.,340.,360.,380.,405.,435.,475.,525.,585.,660.,755.,875.,1025.,1210.,1440.,1730.,2000.,2500.,3200.,4000.])
    h_pred_cl = h_pred.Clone(label)
    h_pred_cl = h_pred_cl.Rebin(len(rebinning)-1,"",rebinning)
    outfile.cd()
    h_pred_cl.Write()


if __name__ == '__main__':

    # load mass distributions
    fpath =OrderedDict()
    fpathPred =OrderedDict()
    massPlotsSignal =OrderedDict()
    mass = OrderedDict()
    mass_plot = OrderedDict()

    regionSignal='SR3'
    regionBckg=''
    systSignal="Gluino"
    searchRegion=regionSignal
    
    if(regionSignal=='SR1'):
        regionBckg='90ias100'
    if(regionSignal=='SR2'):
        regionBckg='99ias100'
    if(regionSignal=='SR3'):
        regionBckg='999ias100'

    os.system('mkdir -p yieldDir')

    text_file_tex_2017 = open('yieldDir/yield'+regionSignal+'_2017.tex', "w")
    text_file_tex_2017.write('\n \\documentclass{article}')
    text_file_tex_2017.write('\n \\begin{document}')
    text_file_tex_2017.write('\n \\begin{center}')
    text_file_tex_2017.write('\n \\begin{tabular}{ |l|c|c|c| } ')
    text_file_tex_2017.write('\n \hline')
    text_file_tex_2017.write('\n Yield '+regionSignal+' & Pred. & Obs. & Signal \\\\')
    text_file_tex_2017.write('\n \hline')
    text_file_tex_2017.write('\n \hline')

    text_file_tex_2018 = open('yieldDir/yield'+regionSignal+'_2018.tex', "w")
    text_file_tex_2018.write('\n \\documentclass{article}')
    text_file_tex_2018.write('\n \\begin{document}')
    text_file_tex_2018.write('\n \\begin{center}')
    text_file_tex_2018.write('\n \\begin{tabular}{ |l|c|c|c| } ')
    text_file_tex_2018.write('\n \hline')
    text_file_tex_2018.write('\n Yield '+regionSignal+' & Pred. & Obs. & Signal \\\\')
    text_file_tex_2018.write('\n \hline')
    text_file_tex_2018.write('\n \hline')

    # load root file
    pathSignal = "/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/"
    pathPred = "/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/"

    yearSignal='2018'
    codeVersionSignal='77p1_Signal'

    ofileBase = rt.TFile("base.root","RECREATE")

    #fpath['Gluino500_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPgluino_M-500_CodeV'+codeVersionSignal+'_v1.root'
    #fpath['Gluino800_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPgluino_M-800_CodeV'+codeVersionSignal+'_v1.root'
    #fpath['Gluino1000_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPgluino_M-1000_CodeV'+codeVersionSignal+'_v1.root'
    #fpath['Gluino1400_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPgluino_M-1400_CodeV'+codeVersionSignal+'_v1.root'
    #fpath['Gluino1600_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPgluino_M-1600_CodeV'+codeVersionSignal+'_v1.root'
    #fpath['Gluino1800_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPgluino_M-1800_CodeV'+codeVersionSignal+'_v1.root'
    #fpath['Gluino2000_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPgluino_M-2000_CodeV'+codeVersionSignal+'_v1.root'
    #fpath['Gluino2200_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPgluino_M-2200_CodeV'+codeVersionSignal+'_v1.root'
    #fpath['Gluino2400_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPgluino_M-2400_CodeV'+codeVersionSignal+'_v1.root'
    #fpath['Gluino2600_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPgluino_M-2600_CodeV'+codeVersionSignal+'_v1.root'
    #fpath['Stop800_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPstop_M-800_CodeV'+codeVersionSignal+'_v1.root'
    #fpath['Stop1000_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPstop_M-1000_CodeV'+codeVersionSignal+'_v1.root'
    #fpath['Stop1200_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPstop_M-1200_CodeV'+codeVersionSignal+'_v1.root'
    #fpath['Stop1400_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPstop_M-1400_CodeV'+codeVersionSignal+'_v1.root'
    #fpath['Stop1600_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPstop_M-1600_CodeV'+codeVersionSignal+'_v1.root'
    #fpath['Stop1800_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPstop_M-1800_CodeV'+codeVersionSignal+'_v1.root'
    #fpath['Stop2000_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPstop_M-2000_CodeV'+codeVersionSignal+'_v1.root'
    #fpath['Stop2200_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPstop_M-2200_CodeV'+codeVersionSignal+'_v1.root'
    #fpath['Stop2400_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPstop_M-2400_CodeV'+codeVersionSignal+'_v1.root'
    #fpath['Stop2600_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPstop_M-2600_CodeV'+codeVersionSignal+'_v1.root'
    #fpath['pairStau308_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPpairStau_M-308_CodeV'+codeVersionSignal+'_v1.root'
    #fpath['pairStau432_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPpairStau_M-432_CodeV'+codeVersionSignal+'_v1.root'
    #fpath['pairStau557_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPpairStau_M-557_CodeV'+codeVersionSignal+'_v1.root'
    #fpath['pairStau651_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPpairStau_M-651_CodeV'+codeVersionSignal+'_v1.root'
    #fpath['pairStau745_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPpairStau_M-745_CodeV'+codeVersionSignal+'_v1.root'
    #fpath['pairStau871_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPpairStau_M-871_CodeV'+codeVersionSignal+'_v1.root'
    #fpath['pairStau1029_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPpairStau_M-1029_CodeV'+codeVersionSignal+'_v1.root'
    #fpath['DYcharge1e_500_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge1e_M-500_CodeV80p0_v1.root'
    #fpath['DYcharge1e_800_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge1e_M-800_CodeV80p0_v1.root'
    #fpath['DYcharge1e_1000_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge1e_M-1000_CodeV80p0_v1.root'
    #fpath['DYcharge1e_1400_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge1e_M-1400_CodeV80p0_v1.root'
    #fpath['DYcharge1e_1800_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge1e_M-1800_CodeV80p0_v1.root'
    #fpath['DYcharge1e_2200_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge1e_M-2200_CodeV80p0_v1.root'
    #fpath['DYcharge1e_2600_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge1e_M-2600_CodeV80p0_v1.root'
    #fpath['DYcharge2e_400_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge2e_M-400_CodeV80p0_v1.root'
    #fpath['DYcharge2e_500_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge2e_M-500_CodeV80p0_v1.root'
    #fpath['DYcharge2e_800_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge2e_M-800_CodeV80p0_v1.root'
    #fpath['DYcharge2e_1000_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge2e_M-1000_CodeV80p0_v1.root'
    #fpath['DYcharge2e_1400_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge2e_M-1400_CodeV80p0_v1.root'
    #fpath['DYcharge2e_1800_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge2e_M-1800_CodeV80p0_v1.root'
    #fpath['DYcharge2e_2200_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge2e_M-2200_CodeV80p0_v1.root'
    #fpath['DYcharge2e_2600_2018'] = pathSignal+'crab_Analysis_'+yearSignal+'_HSCPtauPrimeCharge2e_M-2600_CodeV80p0_v1.root'

    fpath['Gluino500_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPgluino/Merged_HSCPgluino/HSCPgluino_M-500_merged.root'
    fpath['Gluino800_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPgluino/Merged_HSCPgluino/HSCPgluino_M-800_merged.root'
    fpath['Gluino1000_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPgluino/Merged_HSCPgluino/HSCPgluino_M-1000_merged.root'
    fpath['Gluino1400_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPgluino/Merged_HSCPgluino/HSCPgluino_M-1400_merged.root'
    fpath['Gluino1600_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPgluino/Merged_HSCPgluino/HSCPgluino_M-1600_merged.root'
    fpath['Gluino1800_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPgluino/Merged_HSCPgluino/HSCPgluino_M-1800_merged.root'
    fpath['Gluino2000_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPgluino/Merged_HSCPgluino/HSCPgluino_M-2000_merged.root'
    fpath['Gluino2200_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPgluino/Merged_HSCPgluino/HSCPgluino_M-2200_merged.root'
    fpath['Gluino2400_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPgluino/Merged_HSCPgluino/HSCPgluino_M-2400_merged.root'
    fpath['Gluino2600_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPgluino/Merged_HSCPgluino/HSCPgluino_M-2600_merged.root'
    fpath['Stop400_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPstop/Merged_HSCPstop/HSCPstop_M-400_merged.root'
    fpath['Stop500_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPstop/Merged_HSCPstop/HSCPstop_M-500_merged.root'
    fpath['Stop800_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPstop/Merged_HSCPstop/HSCPstop_M-800_merged.root'
    fpath['Stop1000_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPstop/Merged_HSCPstop/HSCPstop_M-1000_merged.root'
    fpath['Stop1200_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPstop/Merged_HSCPstop/HSCPstop_M-1200_merged.root'
    fpath['Stop1400_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPstop/Merged_HSCPstop/HSCPstop_M-1400_merged.root'
    fpath['Stop1600_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPstop/Merged_HSCPstop/HSCPstop_M-1600_merged.root'
    fpath['Stop1800_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPstop/Merged_HSCPstop/HSCPstop_M-1800_merged.root'
    fpath['Stop2000_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPstop/Merged_HSCPstop/HSCPstop_M-2000_merged.root'
    fpath['Stop2200_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPstop/Merged_HSCPstop/HSCPstop_M-2200_merged.root'
    fpath['Stop2400_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPstop/Merged_HSCPstop/HSCPstop_M-2400_merged.root'
    fpath['Stop2600_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPstop/Merged_HSCPstop/HSCPstop_M-2600_merged.root'
    fpath['pairStau308_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPpairStau/Merged_HSCPpairStau/HSCPpairStau_M-308_merged.root'
    fpath['pairStau432_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPpairStau/Merged_HSCPpairStau/HSCPpairStau_M-432_merged.root'
    fpath['pairStau557_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPpairStau/Merged_HSCPpairStau/HSCPpairStau_M-557_merged.root'
    fpath['pairStau651_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPpairStau/Merged_HSCPpairStau/HSCPpairStau_M-651_merged.root'
    fpath['pairStau745_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPpairStau/Merged_HSCPpairStau/HSCPpairStau_M-745_merged.root'
    fpath['pairStau871_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPpairStau/Merged_HSCPpairStau/HSCPpairStau_M-871_merged.root'
    fpath['pairStau1029_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPpairStau/Merged_HSCPpairStau/HSCPpairStau_M-1029_merged.root'
    fpath['DYcharge1e_400_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPtauPrimeCharge1e/Merged_HSCPtauPrimeCharge1e/HSCPtauPrimeCharge1e_M-400_merged.root'
    fpath['DYcharge1e_500_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPtauPrimeCharge1e/Merged_HSCPtauPrimeCharge1e/HSCPtauPrimeCharge1e_M-500_merged.root'
    fpath['DYcharge1e_800_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPtauPrimeCharge1e/Merged_HSCPtauPrimeCharge1e/HSCPtauPrimeCharge1e_M-800_merged.root'
    fpath['DYcharge1e_1000_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPtauPrimeCharge1e/Merged_HSCPtauPrimeCharge1e/HSCPtauPrimeCharge1e_M-1000_merged.root'
    fpath['DYcharge1e_1400_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPtauPrimeCharge1e/Merged_HSCPtauPrimeCharge1e/HSCPtauPrimeCharge1e_M-1400_merged.root'
    fpath['DYcharge1e_1800_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPtauPrimeCharge1e/Merged_HSCPtauPrimeCharge1e/HSCPtauPrimeCharge1e_M-1800_merged.root'
    fpath['DYcharge1e_2200_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPtauPrimeCharge1e/Merged_HSCPtauPrimeCharge1e/HSCPtauPrimeCharge1e_M-2200_merged.root'
    fpath['DYcharge1e_2600_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPtauPrimeCharge1e/Merged_HSCPtauPrimeCharge1e/HSCPtauPrimeCharge1e_M-2600_merged.root'
    fpath['DYcharge2e_400_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPtauPrimeCharge2e/Merged_HSCPtauPrimeCharge2e/HSCPtauPrimeCharge2e_M-400_merged.root'
    fpath['DYcharge2e_500_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPtauPrimeCharge2e/Merged_HSCPtauPrimeCharge2e/HSCPtauPrimeCharge2e_M-500_merged.root'
    #fpath['DYcharge2e_800_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPtauPrimeCharge2e/Merged_HSCPtauPrimeCharge2e/HSCPtauPrimeCharge2e_M-800_merged.root'
    fpath['DYcharge2e_1000_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPtauPrimeCharge2e/Merged_HSCPtauPrimeCharge2e/HSCPtauPrimeCharge2e_M-1000_merged.root'
    fpath['DYcharge2e_1400_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPtauPrimeCharge2e/Merged_HSCPtauPrimeCharge2e/HSCPtauPrimeCharge2e_M-1400_merged.root'
    fpath['DYcharge2e_1800_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPtauPrimeCharge2e/Merged_HSCPtauPrimeCharge2e/HSCPtauPrimeCharge2e_M-1800_merged.root'
    fpath['DYcharge2e_2200_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPtauPrimeCharge2e/Merged_HSCPtauPrimeCharge2e/HSCPtauPrimeCharge2e_M-2200_merged.root'
    fpath['DYcharge2e_2600_2018'] = '/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/Unblinded_hists_dylan/HSCPtauPrimeCharge2e/Merged_HSCPtauPrimeCharge2e/HSCPtauPrimeCharge2e_M-2600_merged.root'


    idirSignal='HSCParticleAnalyzer/BaseName/'
    
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
    nPE='200'

    codeVersion='UnB_v1_v1'
    endLabel='_UnB_v3_Data_v2'

    #codeVersion='73p3'
    #endLabel='_19april_SR3'

    fpathPred['obs_2017'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_nominal'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_etaup'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta2_rebinIh4_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_etadown'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta8_rebinIh4_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_ihup'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh2_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_ihdown'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh8_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_momup'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP1_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_momdown'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP4_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_corrih'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_corrTemplateIh_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_corrmom'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_corrTemplateP_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_fitihup'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_fitIhUp_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_fitihdown'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_fitIhDown_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_fitmomup'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_fitPUp_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_fitmomdown'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_fitPDown_nPE'+nPE+endLabel+'.root')

    year='2018'

    fpathPred['obs_2018'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_nominal'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_etaup'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta2_rebinIh4_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_etadown'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta8_rebinIh4_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_ihup'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh2_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_ihdown'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh8_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_momup'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP1_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_momdown'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP4_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_corrih'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_corrTemplateIh_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_corrmom'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_corrTemplateP_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_fitihup'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_fitIhUp_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_fitihdown'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_fitIhDown_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_fitmomdown'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_fitPDown_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_fitmomup'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_fitPUp_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_fitmomup2'] = rt.TFile(pathPred+'crab_Analysis_SingleMuon_Run'+year+'_CodeV'+codeVersion+'_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_fitPUp_nPE'+nPE+endLabel+'.root')
    
    ###--------------
    # Raph production UnB_V4

    year='2017'
    nPE='200'
    endLabel='_UnB_v4_SR3'

    fpathPred['obs_2017'] = rt.TFile(pathPred+'SingleMuon_'+year+'_Raph_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_nominal'] = rt.TFile(pathPred+'SingleMuon_'+year+'_Raph_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_etaup'] = rt.TFile(pathPred+'SingleMuon_'+year+'_Raph_cutIndex3_rebinEta2_rebinIh4_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_etadown'] = rt.TFile(pathPred+'SingleMuon_'+year+'_Raph_cutIndex3_rebinEta8_rebinIh4_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_ihup'] = rt.TFile(pathPred+'SingleMuon_'+year+'_Raph_cutIndex3_rebinEta4_rebinIh2_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_ihdown'] = rt.TFile(pathPred+'SingleMuon_'+year+'_Raph_cutIndex3_rebinEta4_rebinIh8_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_momup'] = rt.TFile(pathPred+'SingleMuon_'+year+'_Raph_cutIndex3_rebinEta4_rebinIh4_rebinP1_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_momdown'] = rt.TFile(pathPred+'SingleMuon_'+year+'_Raph_cutIndex3_rebinEta4_rebinIh4_rebinP4_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_corrih'] = rt.TFile(pathPred+'SingleMuon_'+year+'_Raph_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_corrTemplateIh_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_corrmom'] = rt.TFile(pathPred+'SingleMuon_'+year+'_Raph_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_corrTemplateP_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_fitihup'] = rt.TFile(pathPred+'SingleMuon_'+year+'_Raph_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_fitIhUp_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_fitihdown'] = rt.TFile(pathPred+'SingleMuon_'+year+'_Raph_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_fitIhDown_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_fitmomup'] = rt.TFile(pathPred+'SingleMuon_'+year+'_Raph_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_fitPUp_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2017_fitmomdown'] = rt.TFile(pathPred+'SingleMuon_'+year+'_Raph_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_fitPDown_nPE'+nPE+endLabel+'.root')

    year='2018'

    fpathPred['obs_2018'] = rt.TFile(pathPred+'SingleMuon_'+year+'_Raph_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_nominal'] = rt.TFile(pathPred+'SingleMuon_'+year+'_Raph_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_etaup'] = rt.TFile(pathPred+'SingleMuon_'+year+'_Raph_cutIndex3_rebinEta2_rebinIh4_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_etadown'] = rt.TFile(pathPred+'SingleMuon_'+year+'_Raph_cutIndex3_rebinEta8_rebinIh4_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_ihup'] = rt.TFile(pathPred+'SingleMuon_'+year+'_Raph_cutIndex3_rebinEta4_rebinIh2_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_ihdown'] = rt.TFile(pathPred+'SingleMuon_'+year+'_Raph_cutIndex3_rebinEta4_rebinIh8_rebinP2_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_momup'] = rt.TFile(pathPred+'SingleMuon_'+year+'_Raph_cutIndex3_rebinEta4_rebinIh4_rebinP1_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_momdown'] = rt.TFile(pathPred+'SingleMuon_'+year+'_Raph_cutIndex3_rebinEta4_rebinIh4_rebinP4_rebinMass1_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_corrih'] = rt.TFile(pathPred+'SingleMuon_'+year+'_Raph_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_corrTemplateIh_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_corrmom'] = rt.TFile(pathPred+'SingleMuon_'+year+'_Raph_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_corrTemplateP_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_fitihup'] = rt.TFile(pathPred+'SingleMuon_'+year+'_Raph_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_fitIhUp_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_fitihdown'] = rt.TFile(pathPred+'SingleMuon_'+year+'_Raph_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_fitIhDown_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_fitmomdown'] = rt.TFile(pathPred+'SingleMuon_'+year+'_Raph_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_fitPDown_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_fitmomup'] = rt.TFile(pathPred+'SingleMuon_'+year+'_Raph_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_fitPUp_nPE'+nPE+endLabel+'.root')
    fpathPred['pred_2018_fitmomup2'] = rt.TFile(pathPred+'SingleMuon_'+year+'_Raph_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_fitPUp_nPE'+nPE+endLabel+'.root')


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
        #p1=0.00131767
        #p0=0.939438
        ### UNBLINDING RAPH
        p1=0.000771218
        p0=0.979842
        
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
    mass_plot['pred_2017_fitihup'] = fpathPred['pred_2017_fitihup'].Get("mass_predBC_"+regionBckg)
    mass_plot['pred_2017_fitihdown'] = fpathPred['pred_2017_fitihdown'].Get("mass_predBC_"+regionBckg)
    mass_plot['pred_2017_fitmomup'] = fpathPred['pred_2017_fitmomup'].Get("mass_predBC_"+regionBckg)
    mass_plot['pred_2017_fitmomdown'] = fpathPred['pred_2017_fitmomdown'].Get("mass_predBC_"+regionBckg)
    
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
    mass_plot['pred_2017_fitihup'] = BiasCorrection(mass_plot['pred_2017_fitihup'],p1,p0)
    mass_plot['pred_2017_fitihdown'] = BiasCorrection(mass_plot['pred_2017_fitihdown'],p1,p0)
    mass_plot['pred_2017_fitmomup'] = BiasCorrection(mass_plot['pred_2017_fitmomup'],p1,p0)
    mass_plot['pred_2017_fitmomdown'] = BiasCorrection(mass_plot['pred_2017_fitmomdown'],p1,p0)
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
        #p1=0.00169109
        #p0=0.901472
        ### UNBLINDING V4 RAPH
        p1=0.00138625
        p0=0.94072

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
    mass_plot['pred_2018_fitihup'] = fpathPred['pred_2018_fitihup'].Get("mass_predBC_"+regionBckg)
    mass_plot['pred_2018_fitihdown'] = fpathPred['pred_2018_fitihdown'].Get("mass_predBC_"+regionBckg)
    mass_plot['pred_2018_fitmomup'] = fpathPred['pred_2018_fitmomup'].Get("mass_predBC_"+regionBckg)
    mass_plot['pred_2018_fitmomdown'] = fpathPred['pred_2018_fitmomdown'].Get("mass_predBC_"+regionBckg)
    mass_plot['pred_2018_fitmomdown'] = fpathPred['pred_2018_fitmomdown'].Get("mass_predBC_"+regionBckg)
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
    mass_plot['pred_2018_fitihup'] = BiasCorrection(mass_plot['pred_2018_fitihup'],p1,p0)
    mass_plot['pred_2018_fitihdown'] = BiasCorrection(mass_plot['pred_2018_fitihdown'],p1,p0)
    mass_plot['pred_2018_fitmomup'] = BiasCorrection(mass_plot['pred_2018_fitmomup'],p1,p0)
    mass_plot['pred_2018_fitmomdown'] = BiasCorrection(mass_plot['pred_2018_fitmomdown'],p1,p0)
    mass_plot['pred_2018_corrbias'] = BiasCorrection(mass_plot['pred_2018_nominal'],p1,p0)

    ofileBase.cd()
    
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
        'Stop400_2018': '$\\tilde{t}$ (M=400 GeV)',
        'Stop500_2018': '$\\tilde{t}$ (M=500 GeV)',
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
       'Stop400_2018': 400,
       'Stop500_2018': 500,
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
    
    outDataCardsDir = "datacards_"+searchRegion+"_test_UnB_v4_Raph_withGoodSignals/"
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

    outMassShape = "mass_shape_analysis_dir"
    os.system("mkdir -p {0}".format(outMassShape))

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

        print 'obs2017: ',obs_2017, ' obs2018: ', obs_2018

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
        bkg_2017_fitIh_up = integralHisto(mass_plot['pred_2017_fitihup'],xmin,xmax)
        bkg_2017_fitIh_down =integralHisto(mass_plot['pred_2017_fitihdown'],xmin,xmax)
        bkg_2017_fitMom_up =integralHisto(mass_plot['pred_2017_fitmomup'],xmin,xmax)
        bkg_2017_fitMom_down =integralHisto(mass_plot['pred_2017_fitmomdown'],xmin,xmax)
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
        bkg_2018_fitIh_up = integralHisto(mass_plot['pred_2018_fitihup'],xmin,xmax)
        bkg_2018_fitIh_down =integralHisto(mass_plot['pred_2018_fitihdown'],xmin,xmax)
        bkg_2018_fitMom_up =integralHisto(mass_plot['pred_2018_fitmomup'],xmin,xmax)
        bkg_2018_fitMom_down =integralHisto(mass_plot['pred_2018_fitmomdown'],xmin,xmax)
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

        ofileSignal_name = signal+".root"
        ofileSignal = rt.TFile(ofileSignal_name,"RECREATE")

        #for systOnSignal in massPlotsSignal.keys():
            #saveSignalShape(ofileSignal, ifile.Get(massPlotsSignal[systOnSignal]), systOnSignal)
        

        #saveSigShape(ofileSignal, ifile.Get(massPlotsSignal['Nominal']), "Nominal")
        #saveSigShape(ofileSignal, ifile.Get(massPlotsSignal['PU_up']), "PUUp")
        #saveSigShape(ofileSignal, ifile.Get(massPlotsSignal['PU_down']), "PUDown")
        #saveSigShape(ofileSignal, ifile.Get(massPlotsSignal['Fpix_up']), "FpixUp")
        #saveSigShape(ofileSignal, ifile.Get(massPlotsSignal['Fpix_down']), "FpixDown")
        #saveSigShape(ofileSignal, ifile.Get(massPlotsSignal['Gstrip_up']), "GstripUp")
        #saveSigShape(ofileSignal, ifile.Get(massPlotsSignal['Gstrip_down']), "GstripDown")
        #saveSigShape(ofileSignal, ifile.Get(massPlotsSignal['Pt_up']), "PtUp")
        #saveSigShape(ofileSignal, ifile.Get(massPlotsSignal['Pt_down']), "PtDown")
        #saveSigShape(ofileSignal, ifile.Get(massPlotsSignal['Trigger_up']), "TriggerUp")
        #saveSigShape(ofileSignal, ifile.Get(massPlotsSignal['Trigger_down']), "TriggerDown")
        #saveSigShape(ofileSignal, ifile.Get(massPlotsSignal['K_up']), "KUp")
        #saveSigShape(ofileSignal, ifile.Get(massPlotsSignal['K_down']), "KDown")
        #saveSigShape(ofileSignal, ifile.Get(massPlotsSignal['C_up']), "CUp")
        #saveSigShape(ofileSignal, ifile.Get(massPlotsSignal['C_down']), "CDown")
#
        #savePredShape(ofileSignal, mass_plot['obs_2017'], "data_obs_Ch2017")
        #savePredShape(ofileSignal, mass_plot['obs_2018'], "data_obs_Ch2018")
#
        #savePredShape(ofileSignal, mass_plot['pred_2017_nominal'], "background_Ch2017")
        #savePredShape(ofileSignal, mass_plot['pred_2017_stat'], "background_Ch2017_statUp")
        #savePredShape(ofileSignal, mass_plot['pred_2017_stat'], "background_Ch2017_statDown")
        #savePredShape(ofileSignal, mass_plot['pred_2017_etaup'], "background_Ch2017_eta_up")
        #savePredShape(ofileSignal, mass_plot['pred_2017_etadown'], "background_Ch2017_eta_down")
        ##savePredShape(ofileSignal, mass_plot['pred_2017_etadown'], "background_Ch2017_etaUp")
        ##savePredShape(ofileSignal, mass_plot['pred_2017_etadown'], "background_Ch2017_etaDown")
        #savePredShape(ofileSignal, mass_plot['pred_2017_ihup'], "background_Ch2017_ih_up")
        #savePredShape(ofileSignal, mass_plot['pred_2017_ihdown'], "background_Ch2017_ih_down")
        ##savePredShape(ofileSignal, mass_plot['pred_2017_ihdown'], "background_Ch2017_ihUp")
        ##savePredShape(ofileSignal, mass_plot['pred_2017_ihdown'], "background_Ch2017_ihDown")
        #savePredShape(ofileSignal, mass_plot['pred_2017_momup'], "background_Ch2017_mom_up")
        #savePredShape(ofileSignal, mass_plot['pred_2017_momdown'], "background_Ch2017_mom_down")
        ##savePredShape(ofileSignal, mass_plot['pred_2017_momdown'], "background_Ch2017_momUp")
        ##savePredShape(ofileSignal, mass_plot['pred_2017_momdown'], "background_Ch2017_momDown")
        #savePredShape(ofileSignal, mass_plot['pred_2017_corrih'], "background_Ch2017_corrihUp")
        #savePredShape(ofileSignal, mass_plot['pred_2017_corrih'], "background_Ch2017_corrihDown")
        #savePredShape(ofileSignal, mass_plot['pred_2017_corrmom'], "background_Ch2017_corrmomUp")
        #savePredShape(ofileSignal, mass_plot['pred_2017_corrmom'], "background_Ch2017_corrmomDown")
        #savePredShape(ofileSignal, mass_plot['pred_2017_fitihup'], "background_Ch2017_fitihUp")
        #savePredShape(ofileSignal, mass_plot['pred_2017_fitihdown'], "background_Ch2017_fitihDown")
        #savePredShape(ofileSignal, mass_plot['pred_2017_fitmomup'], "background_Ch2017_fitmomUp")
        #savePredShape(ofileSignal, mass_plot['pred_2017_fitmomdown'], "background_Ch2017_fitmomDown")
        #savePredShape(ofileSignal, mass_plot['pred_2017_corrbias'], "background_Ch2017_biascorrectionUp")
        #savePredShape(ofileSignal, mass_plot['pred_2017_corrbias'], "background_Ch2017_biascorrectionDown")
        #
        #savePredShape(ofileSignal, mass_plot['pred_2018_nominal'], "background_Ch2018")
        #savePredShape(ofileSignal, mass_plot['pred_2018_stat'], "background_Ch2018_statUp")
        #savePredShape(ofileSignal, mass_plot['pred_2018_stat'], "background_Ch2018_statDown")
        #savePredShape(ofileSignal, mass_plot['pred_2018_etaup'], "background_Ch2018_eta_up")
        #savePredShape(ofileSignal, mass_plot['pred_2018_etadown'], "background_Ch2018_eta_down")
        ##savePredShape(ofileSignal, mass_plot['pred_2018_etadown'], "background_Ch2018_etaUp")
        ##savePredShape(ofileSignal, mass_plot['pred_2018_etadown'], "background_Ch2018_etaDown")
        #savePredShape(ofileSignal, mass_plot['pred_2018_ihup'], "background_Ch2018_ih_up")
        #savePredShape(ofileSignal, mass_plot['pred_2018_ihdown'], "background_Ch2018_ih_down")
        ##savePredShape(ofileSignal, mass_plot['pred_2018_ihdown'], "background_Ch2018_ihUp")
        ##savePredShape(ofileSignal, mass_plot['pred_2018_ihdown'], "background_Ch2018_ihDown")
        #savePredShape(ofileSignal, mass_plot['pred_2018_momup'], "background_Ch2018_mom_up")
        #savePredShape(ofileSignal, mass_plot['pred_2018_momdown'], "background_Ch2018_mom_down")
        ##savePredShape(ofileSignal, mass_plot['pred_2018_momdown'], "background_Ch2018_momUp")
        ##savePredShape(ofileSignal, mass_plot['pred_2018_momdown'], "background_Ch2018_momDown")
        #savePredShape(ofileSignal, mass_plot['pred_2018_corrih'], "background_Ch2018_corrihUp")
        #savePredShape(ofileSignal, mass_plot['pred_2018_corrih'], "background_Ch2018_corrihDown")
        #savePredShape(ofileSignal, mass_plot['pred_2018_corrmom'], "background_Ch2018_corrmomUp")
        #savePredShape(ofileSignal, mass_plot['pred_2018_corrmom'], "background_Ch2018_corrmomDown")
        #savePredShape(ofileSignal, mass_plot['pred_2018_fitihup'], "background_Ch2018_fitihUp")
        #savePredShape(ofileSignal, mass_plot['pred_2018_fitihdown'], "background_Ch2018_fitihDown")
        #savePredShape(ofileSignal, mass_plot['pred_2018_fitmomup'], "background_Ch2018_fitmomUp")
        #savePredShape(ofileSignal, mass_plot['pred_2018_fitmomdown'], "background_Ch2018_fitmomDown")
        #savePredShape(ofileSignal, mass_plot['pred_2018_corrbias'], "background_Ch2018_biascorrectionUp")
        #savePredShape(ofileSignal, mass_plot['pred_2018_corrbias'], "background_Ch2018_biascorrectionDown")
        #os.system("mv {0} {1}/.".format(ofileSignal_name,outMassShape))

        yield_total=0
        yield_minus2sigma_minus1sigma=0
        yield_minus1sigma_mu=0
        yield_mu_plus1sigma=0
        yield_plus1sigma_plus2sigma=0
        yield_minus1sigma_plus2sigma=0

        #print (int(mean-stddev), int(mean+2*stddev))
        #print ("integral 2017: ", mass_plot['obs_2017'].Integral(mass_plot['obs_2017'].FindBin(mean-stddev),mass_plot['obs_2017'].FindBin(mean+2*stddev)))
        #print ("integral 2018: ", mass_plot['obs_2018'].Integral(mass_plot['obs_2018'].FindBin(mean-stddev),mass_plot['obs_2018'].FindBin(mean+2*stddev)))

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
            print (signal_Trigger_up,signal_Trigger_down)
            print pushSyst(signal_Trigger_up,signal_Trigger_down)
            yTrigger.append(pushSyst(signal_Trigger_up,signal_Trigger_down))
            yK.append(pushSyst(signal_K_up,signal_K_down))
            yC.append(pushSyst(signal_C_up,signal_C_down))
            print np.sqrt(pow(pushSyst(signal_pu_up,signal_pu_down),2)+pow(pushSyst(signal_Fpix_up,signal_Fpix_down),2)+pow(pushSyst(signal_Gstrip_up,signal_Gstrip_down),2)+pow(pushSyst(signal_Pt_up,signal_Pt_down),2)+pow(pushSyst(signal_Trigger_up,signal_Trigger_down),2)+pow(pushSyst(signal_K_up,signal_K_down),2)+pow(pushSyst(signal_C_up,signal_C_down),2))
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
            bkg_2017_fitIh_up/=bkg_2017_nominal
            bkg_2017_fitIh_down/=bkg_2017_nominal
            bkg_2017_fitMom_up/=bkg_2017_nominal
            bkg_2017_fitMom_down/=bkg_2017_nominal
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
            bkg_2018_fitIh_up/=bkg_2018_nominal
            bkg_2018_fitIh_down/=bkg_2018_nominal
            bkg_2018_fitMom_up/=bkg_2018_nominal
            bkg_2018_fitMom_down/=bkg_2018_nominal
            bkg_2018_correctionBias/=bkg_2018_nominal

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

        sig_unc_correlated = {
            'sig_pu' : [signal_pu_down, signal_pu_up, signal_pu_down, signal_pu_up],
            'sig_Trigger' : [signal_Trigger_down, signal_Trigger_up, signal_Trigger_down, signal_Trigger_up],
        }

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

        bkg_unc_correlated = {
            'bkg_etaBinning' : [bkg_2017_etaBinning_down, bkg_2017_etaBinning_up, bkg_2018_etaBinning_down, bkg_2018_etaBinning_up],
            'bkg_ihBinning' : [bkg_2017_ihBinning_down, bkg_2017_ihBinning_up, bkg_2018_ihBinning_down, bkg_2018_ihBinning_up],
            'bkg_momentumBinning' : [bkg_2017_momBinning_down, bkg_2017_momBinning_up, bkg_2018_momBinning_down, bkg_2018_momBinning_up],
            'bkg_correctionTemplateIh' : [bkg_2017_corrTemplateIh_up, bkg_2018_corrTemplateIh_up],
            'bkg_correctionTemplateMomentum' : [bkg_2017_corrTemplateMom_up, bkg_2018_corrTemplateMom_up],
            'bkg_correctionBias' : [bkg_2017_correctionBias, bkg_2018_correctionBias],
        }

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

        signal_yield_norm=signal_yield*41.5/101
        signal_yield_norm=signal_yield*59.7/101

        ### fill the datacard
        make_datacard_hscp_combining2017and2018(outDataCardsDir, signal, signal_yield*41.5/101., signal_yield*59.7/101., bkg_2017_nominal, bkg_2018_nominal, obs_2017, obs_2018, sig_2017_unc, sig_2018_unc, sig_unc_correlated, bkg_2017_unc, bkg_2018_unc, bkg_unc_correlated)
        
        bkg=bkg_2017_nominal
        bckg_up=totalUncertainy(tot_unc_bkg_2017,1)
        bckg_down=totalUncertainy(tot_unc_bkg_2017,0)
        signal_up=totalUncertainy(tot_unc_sig_2017,1)
        signal_down=totalUncertainy(tot_unc_sig_2017,0)
        signal_yield_norm=signal_yield*41.5/101.
        
        makeYieldFile(text_file_tex_2017,name[signal],bkg,bckg_up,bckg_down,obs_2017,signal_yield_norm,signal_up,signal_down)
        
        bkg=bkg_2018_nominal
        bckg_up=totalUncertainy(tot_unc_bkg_2018,1)
        bckg_down=totalUncertainy(tot_unc_bkg_2018,0)
        signal_up=totalUncertainy(tot_unc_sig_2018,1)
        signal_down=totalUncertainy(tot_unc_sig_2018,0)
        signal_yield_norm=signal_yield*59.7/101.
        
        makeYieldFile(text_file_tex_2018,name[signal],bkg,bckg_up,bckg_down,obs_2018,signal_yield_norm,signal_up,signal_down)

    grPU = create_TGraph(x, yPU, axis_title=['target mass [GeV]', 'Systematics uncertainties [%]'])
    grFpix = create_TGraph(x, yFpix, axis_title=['target mass [GeV]', 'Systematics uncertainties [%]'])
    grGstrip = create_TGraph(x, yGstrip, axis_title=['target mass [GeV]', 'Systematics uncertainties [%]'])
    grPt = create_TGraph(x, yPt, axis_title=['target mass [GeV]', 'Systematics uncertainties [%]'])
    grTrigger = create_TGraph(x, yTrigger, axis_title=['target mass [GeV]', 'Systematics uncertainties [%]'])
    grK = create_TGraph(x, yK, axis_title=['target mass [GeV]', 'Systematics uncertainties [%]'])
    grC = create_TGraph(x, yC, axis_title=['target mass [GeV]', 'Systematics uncertainties [%]'])
    grTotal = create_TGraph(x, yTotal, axis_title=['target mass [GeV]', 'Systematics uncertainties [%]'])


    text_file_tex_2017.write('\n \\end{tabular}')
    text_file_tex_2017.write('\n \\end{center}')
    text_file_tex_2017.write('\n \\end{document}')
    text_file_tex_2017.close()

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
    CMS_lumi.lumi_13TeV  = "101 fb^{-1}"
    
    rt.gStyle.SetOptStat(0)
    c1=rt.TCanvas()
    h2.Draw("col")
    h2.GetXaxis().SetTitleSize(0.04)
    h2.GetYaxis().SetTitleSize(0.04)
    h2.GetXaxis().SetTitleOffset(1.5)
    h2.GetYaxis().SetTitleOffset(2)
    CMS_lumi.CMS_lumi(c1, 4, iPos)
    os.system('mkdir -p mass_window_dir')
    c1.SaveAs("mass_window_dir/h2_massWindow"+regionSignal+".root")
    c1.SaveAs("mass_window_dir/h2_massWindow"+regionSignal+".pdf")

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
    os.system('mkdir -p systTargetMass_dir')
    c2.SaveAs("systTargetMass_dir/systTargetMass_"+systSignal+"_"+regionSignal+".root")
    c2.SaveAs("systTargetMass_dir/systTargetMass_"+systSignal+"_"+regionSignal+".pdf")
