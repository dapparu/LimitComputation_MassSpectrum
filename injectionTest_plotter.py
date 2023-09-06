import ROOT as rt
import csv
import re
import sys
import collections
import os

from collections import OrderedDict
import uproot

import scipy
import awkward
import numpy as np
import time

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i","--inputFile",type=str,help="input file path")
parser.add_argument("-o","--outputFile",type=str,help="output file name")
parser.add_argument("-exps","--expectedsignal",type=float,help="expected signal to inject")
parser.add_argument("-fL","--fitLikelihood",type=int,help="gaussian likelihood fit")
parser.add_argument("-r","--rangeBins",type=int,help="range of the pull distribution")
parser.add_argument("-n","--nbins",type=int,help="number of bins of the pull distribution")
args = parser.parse_args()

print "------Doing signal injection test plot------"

ifile = args.inputFile
ofile = args.outputFile
exps  = args.expectedsignal
print (ifile.split("expected")[0])

likelihoodFit = args.fitLikelihood
if(args.rangeBins>0): rb = args.rangeBins
else: rb = 5
if(args.nbins>0): nb = args.nbins
else: nb=20

rt.gROOT.SetBatch(True)
rt.gStyle.SetOptStat(True)
rt.gStyle.SetOptFit(True)


injectedAmount = exps

injectedName = str(injectedAmount)
toyOutput = rt.TFile.Open(ifile)

tree_fit_sb = toyOutput.Get('tree_fit_sb')

# Final plotting
result_can = rt.TCanvas('sigpull_can','sigpull_can',800,800)
#tree_fit_sb.Draw("(r-{rinj})/(rHiErr*(r<{rinj})+rLoErr*(r>{rinj}))>>sigpull({nbs},-{rbs},{rbs})".format(rinj=injectedAmount,nbs=nb,rbs=rb),"fit_status>=0")
tree_fit_sb.Draw("(r-{rinj})/((rHiErr+rLoErr)/2.)>>sigpull({nbs},-{rbs},{rbs})".format(rinj=injectedAmount,nbs=nb,rbs=rb),"fit_status>=0")
hsigpull = rt.gDirectory.Get('sigpull')
if(likelihoodFit==1):  
    hsigpull.Fit("gaus","L") #log-likelihood fit
else:
    hsigpull.Fit("gaus","")
hsigpull.SetTitle('')
hsigpull.GetXaxis().SetTitle('(r-%s)/rErr'%injectedAmount)
hsigpull.GetYaxis().SetTitle('Yield')
hsigpull.GetYaxis().SetTitleOffset(1.2)
result_can.cd()
hsigpull.Draw('pe')

print("------------------ coverage test ------------------")
coverage_test = hsigpull.Integral(hsigpull.FindBin(-1),hsigpull.FindBin(1)-1)/hsigpull.Integral()
print('coverage: ',coverage_test)

result_can.Print(ofile+'_r%s_pull.pdf'%(str(injectedAmount).replace('.','p')),'pdf')
result_can.Print(ofile+'_r%s_pull.root'%(str(injectedAmount).replace('.','p')),'root')
result_can.Print(ofile+'_r%s_pull.C'%(str(injectedAmount).replace('.','p')),'C')

