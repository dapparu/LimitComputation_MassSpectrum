#This is the script that takes a datacard and an injected signal strenght
#this makes the toys, and does the Fit Diagnostics

import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-d","--datacard",type=str,help="data card to input")
parser.add_argument("-i","--info",type=str,help="extra info to append to name")
parser.add_argument("-exps","--expectedsignal",type=float,help="expected signal to inject")
parser.add_argument("-t","--ntoys",type=int,help="number of toys to generate")
parser.add_argument("-a","--fitAlgo",type=str,help="set robust fit algorithm used by combine")
args = parser.parse_args()

print "------Doing signal injection test------"

dcard = args.datacard
info  = args.info
exps  = args.expectedsignal
ntoys = args.ntoys

algo_name = "Minuit2,Combined"

algo_name = args.fitAlgo

splitnz = dcard.split(".txt")[0]
name = splitnz
print "The name signifier for the output is ",name+"_expectedsignal"+str(exps)+"_"+info+".txt"
nametocombine = name+"_expectedsignal"+str(exps)+"_"+info
#Generate the toys
print "   ---Generating Injected Toys---   "
print "   injecting signal r = ",exps
print "   generating toys: ",ntoys

#The combine incantation
# do NOT use --toysFrequentist --bypassFrequentist options because we want to take into account systematics errors in the toy production

subprocess.call("combine -M GenerateOnly --saveToys -d "+dcard+" -n "+nametocombine+" -t "+str(ntoys)+" --expectSignal "+str(exps),shell=True)
print("combine -M GenerateOnly --saveToys -d "+dcard+" -n "+nametocombine+" -t "+str(ntoys)+" --expectSignal "+str(exps))
print "     ---Doing Fit Diagnositcs---   "

if("inflatedBkg" in dcard):
	subprocess.call("combine -M FitDiagnostics --plots --setRobustFitAlgo="+algo_name+" --robustFit=1 -d "+dcard+" -n "+nametocombine+" -t "+str(ntoys)+" --rMin -5.0 --rMax 5.0 --toysFile=higgsCombine"+nametocombine+".GenerateOnly.mH120.123456.root"+" --freezeParameters bkgRate",shell=True)
	print("combine -M FitDiagnostics -d "+dcard+" -n "+nametocombine+" -t "+str(ntoys)+" --rMin -5.0 --rMax 5.0 --toysFile=higgsCombine"+nametocombine+".GenerateOnly.mH120.123456.root"+" --freezeParameters bkgRate")
else:
	subprocess.call("combine -M FitDiagnostics --plots --setRobustFitAlgo="+algo_name+" --robustFit=1 -d "+dcard+" -n "+nametocombine+" -t "+str(ntoys)+" --rMin -20.0 --rMax 20.0 --toysFile=higgsCombine"+nametocombine+".GenerateOnly.mH120.123456.root",shell=True)
	print("combine -M FitDiagnostics --plots --setRobustFitAlgo="+algo_name+" --robustFit=1 -d "+dcard+" -n "+nametocombine+" -t "+str(ntoys)+" --rMin -20.0 --rMax 20.0 --toysFile=higgsCombine"+nametocombine+".GenerateOnly.mH120.123456.root")
	subprocess.call("combine -M Significance -d "+dcard+" -n "+nametocombine+" -t "+str(ntoys)+" --rMin -20.0 --rMax 20.0 --toysFile=higgsCombine"+nametocombine+".GenerateOnly.mH120.123456.root",shell=True)
	print("combine -M Significance -d "+dcard+" -n "+nametocombine+" -t "+str(ntoys)+" --rMin -20.0 --rMax 20.0 --toysFile=higgsCombine"+nametocombine+".GenerateOnly.mH120.123456.root")
	
