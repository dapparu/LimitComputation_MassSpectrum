#! /usr/bin/env python
# 07/02/23 by oponcet

import ROOT; ROOT.PyConfig.IgnoreCommandLineOptions = True
import os, sys, re
import yaml
import CombineHarvester.CombineTools.ch as ch
from CombineHarvester.CombineTools.ch import CombineHarvester, MassesFromRange, SystMap, BinByBinFactory, CardWriter, SetStandardBinNames, AutoRebin
import CombineHarvester.CombinePdfs.morphing as morphing
from CombineHarvester.CombinePdfs.morphing import BuildRooMorphing
from CombineHarvester.CombinePdfs.morphing import BuildCMSHistFuncFactory

import ROOT
from ROOT import RooWorkspace, TFile, RooRealVar


def harvest(setup,**kwargs):
    """Harvest cards."""

    indir       = kwargs.get('indir',     'input' )
    outdir      = kwargs.get('outdir',    'outdir') 
    verbosity   = kwargs.get('verbosity', '1'     )
    analysis    = kwargs.get('analysis',  'HSCParticleAnalyzer'   )
    era         = kwargs.get('era',       '13TeV' )

    

    analysis = setup["analysis"]
    channel = setup["channel"]
    signals = []
    backgrounds = []

    # Categrorise process if it's singal and background 
    for proc in setup["processes"]:
        if "PostS_SR3_Mass" in  proc :
            signals.append(proc)
        elif not "obs" in proc:
            backgrounds.append(proc)

    print "Signals: %s"%signals
    print "Backgrounds: %s"%backgrounds

    cats = [(1,channel)] #region = bin = dans mon cas mes DM

    print "cats: ", cats

    harvester = CombineHarvester()

    harvester.AnalysisSet(analysis)

    # Change flag causing bug : 
    harvester.SetFlag("workspaces-use-clone", True)

    harvester.AddObservations(['*'], [analysis], [era], [channel])
    harvester.AddProcesses(['*'], [analysis], [era], [channel], backgrounds, cats, False)
    harvester.AddProcesses(['*'], [analysis], [era], [channel], signals, cats, True)
    
    filename = "datacard_massShape.root"

    ifileSignal     = "crab_Analysis_2018_HSCPpairStau_M-557_CodeV77p1_Signal_v1.root"
    ifileBckgObs    = "crab_Analysis_SingleMuon_Run2017_2018_CodeV73p3_v4_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_nPE200_25april_VR_SR_cp.root"

    # Add systematics
    if "systematics" in setup:
        for sys in setup["systematics"]:
            sysDef = setup["systematics"][sys]
            scaleFactor = 1.0  
            if "scaleFactor" in sysDef:
                scaleFactor = sysDef["scaleFactor"]
            harvester.cp().process(sysDef["processes"]).AddSyst(harvester, sysDef["name"] if "name" in sysDef else sys, sysDef["effect"], SystMap()(scaleFactor))

    # EXTRACT SHAPES
    print(">>> extracting shapes...")
    #filename = "%s/input.root"%(indir)
    #print ">>>   file %s"%(filename)
    #harvester.cp().channel([channel]).backgrounds().ExtractShapes(ifileBckgObs, "$PROCESS", "")
    harvester.cp().channel([channel]).backgrounds().ExtractShapes(ifileBckgObs, "$PROCESS", "$PROCESS_$SYSTEMATIC")
    harvester.cp().channel([channel]).signals().ExtractShapes(ifileSignal, "$ANALYSIS/$BIN/$PROCESS", "$ANALYSIS/$BIN/$PROCESS_$SYSTEMATIC")
    #harvester.cp().channel([channel]).signals().ExtractShapes(ifileSignal, "HSCParticleAnalyzer/BaseName/$PROCESS", "HSCParticleAnalyzer/BaseName/$PROCESS_$SYSTEMATIC")
    harvester.cp().channel([channel]).signals().ExtractShapes(ifileSignal, "HSCParticleAnalyzer/BaseName/$PROCESS", "HSCParticleAnalyzer/BaseName/$PROCESS_$SYSTEMATIC")



    # ROOVAR
    workspace = RooWorkspace(analysis,analysis)
    #BuildRooMorphing(workspace, harvester, channel, "PostS_SR3_Mass")

    #workspace.Print()
    workspace.writeToFile("workspace_py.root")

    print(analysis)

    # EXTRACT PDFs
    print(">>> add workspace and extract pdf...")
    harvester.AddWorkspace(workspace, False)  
    
    #harvester.ExtractPdfs(harvester, channel, "$PROCESS", "")  # Extract all processes (signal and bkg are named the same way)    
    #harvester.ExtractData(channel, "data_obs")  # Extract the RooDataHist

    #PRINT
    if verbosity>0:
        print("\n>>> print observation...\n")
        harvester.PrintObs()
        print("\n>>> print processes...\n")
        harvester.PrintProcs()
        print("\n>>> print systematics...\n")
        harvester.PrintSysts()
        print("\n>>> print parameters...\n")
        harvester.PrintParams()
        print "\n"

    

    # WRITER
    print(">>> writing datacards...")
    datacardtxt  = "datacard_simple_%s.txt"%(channel)
    datacardroot = "datacard_simple_%s.root"%(channel)
    print datacardtxt, datacardroot
    harvester.WriteDatacard(datacardtxt,datacardroot)
    writer = CardWriter(datacardtxt,datacardroot)
    writer.SetVerbosity(verbosity)
    writer.SetWildcardMasses([ ])
    #writer.WriteCards("outdir", harvester)


def main(args):

    ## Open and import information from config file here to be publicly accessible in all functions
    print "Using configuration file: %s"%args.config
    with open(args.config, 'r') as file:
        setup = yaml.safe_load(file)
    verbosity = 1 if args.verbose else 0
    verbosity = 1 
    indir = ""
    harvest(setup,indir=indir,verbosity=verbosity, analysis='2018')
    
if __name__ == '__main__':

  from argparse import ArgumentParser
  argv = sys.argv
  description = '''This script makes datacards with CombineHarvester.'''
  parser = ArgumentParser(prog="harvesterDatacards",description=description,epilog="Succes!")
  parser.add_argument('-c', '--config', dest='config', type=str, default='/config/configfile.yml', action='store', help="set config file containing sample & fit setup")
  parser.add_argument('-v', '--verbose', dest='verbose', default=False, action='store_true', help="set verbose")
  args = parser.parse_args()

  main(args)
  print ">>>\n>>> done harvesting\n"
    
