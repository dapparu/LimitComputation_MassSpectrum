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


fpath =OrderedDict()
fpathPred =OrderedDict()
massPlotsSignal =OrderedDict()
mass = OrderedDict()
mass_plot = OrderedDict()

regionSignal='SR3'
regionBckg=''

idirSignal='HSCParticleAnalyzer/BaseName/'
searchRegion=regionSignal
regionBckg=''

pathPred = "/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/"

if(searchRegion=='SR1'):
    regionBckg='90ias100'
if(searchRegion=='SR2'):
    regionBckg='99ias100'
if(searchRegion=='SR3'):
    regionBckg='999ias100'


year='2017'
codeVersion='73p3_v4'
nPE='200'
endLabel='_19april'
if(searchRegion=='SR3'):
    endLabel='_21april_SR3'

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

mass_plot['pred_2017_nominal'] = fpathPred['pred_2017_nominal'].Get("mass_predBC_"+regionBckg)
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
mass_plot['pred_2017_stat'] = statErr(mass_plot['pred_2017_nominal'],"statErr")


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
mass_plot['pred_2018_stat'] = statErr(mass_plot['pred_2018_nominal'],"statErr")

ofile = rt.TFile()

print mass_plot['pred_2018_nominal'].Integral()
