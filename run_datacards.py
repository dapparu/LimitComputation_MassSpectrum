"""
Code to combine the datacard of 3 categories in 2tag and produce limits for 2tag limits
"""
import numpy as np
import os
import sys

from threading import Thread

def readNorm(f_cscCard):
    f = open(f_cscCard,"r")
    norm = float(f.readline().split()[3])
    return norm

def insert(originalfile,string):
    with open(originalfile,'r') as f:
        with open('newfile.txt','w') as f2:
            f2.write(string)
            f2.write(f.read())
    os.rename('newfile.txt',originalfile)

if __name__ == '__main__':
	#input_dir=os.getenv('CMSSW_BASE')+'/src/HiggsAnalysis/HSCPLimit/combine/datacards/mass/v_SR1_opti/'
	#tree_dir=os.getenv('CMSSW_BASE')+'/src/HiggsAnalysis/HSCPLimit/combine/limitTrees/mass/v_SR1_opti/'
    
    #input_dir='datacards_test_SR1_5may_v3'
    #tree_dir='limitTrees_test_SR1_5may_v3'

    #input_dir='datacards_test_SR3_8june_noSignalSystematics'
    #tree_dir='limitTrees_test_SR3_8june_noSignalSystematics'

    input_dir='datacards_test_SR3_rescalePAS2016_MuTriggerEfficiency'
    tree_dir='limitTrees_test_SR3_rescalePAS2016_MuTriggerEfficiency'

    input_dir='datacards_test_SR3_rescalePAS2016'
    tree_dir='limitTrees_test_SR3_rescalePAS2016'

    input_dir='datacards_test_SR3_July11_DYToTauPrime'
    tree_dir='limitTrees_test_SR3_July11_DYToTauPrime'

    input_dir='datacards_UnB_v1_SR3_Aug29'
    tree_dir='limitTrees_UnB_v1_SR3_Aug29'
	
    os.system("mkdir -p {0}/".format(input_dir))
    os.system("mkdir -p {0}/".format(tree_dir))
    '''
    samples = [
        ##'Gluino400',
        'Gluino500',
        'Gluino800',
        'Gluino1000',
        ##'Gluino1200',
        'Gluino1400',
        'Gluino1600',
	    'Gluino1800',
        'Gluino2000',
        'Gluino2200',
        'Gluino2400',
        'Gluino2600',
        ##'Stop200',
        ##'Stop600',
        'Stop800',
        'Stop1000',
        'Stop1200',
        'Stop1400',
        'Stop1600',
	    'Stop1800',
        'Stop2000',
        'Stop2200',
        'Stop2400',
        'Stop2600',
        'pairStau308',
        'pairStau432',
        'pairStau557',
        'pairStau651',
        'pairStau745',
        'pairStau871',
        'pairStau1029',
        ##'pairStau308',
        ##'pairStau494',
        ##'pairStau651',
        ##'pairStau1029',
        #'gmsbStau432',
        #'gmsbStau557',
        #'gmsbStau651',
        #'gmsbStau745',
        #'gmsbStau871',
        #'DYcharge1e_500',
        #'DYcharge1e_800',
        #'DYcharge1e_1000',
        #'DYcharge1e_1400',
        #'DYcharge1e_1800',
        #'DYcharge1e_2200',
        #'DYcharge1e_2600',
        #'DYcharge2e_500',
        #'DYcharge2e_800',
        #'DYcharge2e_1000',
        #'DYcharge2e_1400',
        #'DYcharge2e_1800',
        #'DYcharge2e_2200',
        #'DYcharge2e_2600',

    ]
    '''

    samples = [
        ##'Gluino400',
        'Gluino500_2018',
        'Gluino800_2018',
        'Gluino1000_2018',
        ##'Gluino1200',
        'Gluino1400_2018',
        'Gluino1600_2018',
	    'Gluino1800_2018',
        'Gluino2000_2018',
        'Gluino2200_2018',
        'Gluino2400_2018',
        'Gluino2600_2018',
        ##'Stop200',
        ##'Stop600',
        'Stop800_2018',
        'Stop1000_2018',
        'Stop1200_2018',
        'Stop1400_2018',
        'Stop1600_2018',
	    'Stop1800_2018',
        'Stop2000_2018',
        'Stop2200_2018',
        'Stop2400_2018',
        'Stop2600_2018',
        'pairStau308_2018',
        'pairStau432_2018',
        'pairStau557_2018',
        'pairStau651_2018',
        'pairStau745_2018',
        'pairStau871_2018',
        'pairStau1029_2018',
        ##'pairStau308',
        ##'pairStau494',
        ##'pairStau651',
        ##'pairStau1029',
        #'gmsbStau432',
        #'gmsbStau557',
        #'gmsbStau651',
        #'gmsbStau745',
        #'gmsbStau871',
        'DYcharge1e_500_2018',
        'DYcharge1e_800_2018',
        'DYcharge1e_1000_2018',
        'DYcharge1e_1400_2018',
        'DYcharge1e_1800_2018',
        'DYcharge1e_2200_2018',
        'DYcharge1e_2600_2018',
        'DYcharge2e_500_2018',
        'DYcharge2e_800_2018',
        'DYcharge2e_1000_2018',
        'DYcharge2e_1400_2018',
        'DYcharge2e_1800_2018',
        'DYcharge2e_2200_2018',
        'DYcharge2e_2600_2018',

    ]
    
    #for sample in samples:
    def task(sample):
        name=sample
        print name
        #run_combine = "combine -M AsymptoticLimits --run blind -v 1 -n .{} -d {}/{}_2018.txt --rRelAcc 0.000005 --rAbsAcc 0.000005".format(name, input_dir, name)
        #run_combine = "combine -M AsymptoticLimits --run blind -v 1 -n .{} -d {}/{}.txt --rRelAcc 0.000005 --rAbsAcc 0.000005".format(name, input_dir, name)
        run_combine = "combine -M AsymptoticLimits -n .{} -d {}/{}.txt --rRelAcc 0.000005 --rAbsAcc 0.000005 --rMin -1000.0 --rMax 1000.0".format(name, input_dir, name)
        #run_combine = "combine -M AsymptoticLimits --run blind -v 1 -n .{} -d {}/{}_rescaledPAS2016.txt --rRelAcc 0.000005 --rAbsAcc 0.000005".format(name, input_dir, name)
        print run_combine
        os.system(run_combine)
        #run_combine = "combine -M Significance -n.{} {}/{}_2018.txt -t -1 --expectSignal=1".format(name, input_dir, name)
        #run_combine = "combine -M Significance -n.{} {}/{}.txt -t -1 --expectSignal=1".format(name, input_dir, name)
        run_combine = "combine -M Significance -n.{} {}/{}.txt".format(name, input_dir, name)
        #run_combine = "combine -M Significance -n.{} {}/{}_rescaledPAS2016.txt -t -1 --expectSignal=1".format(name, input_dir, name)
        print run_combine
        os.system(run_combine)
        os.system("mv higgsCombine.{0}.AsymptoticLimits.mH120.root {1}/".format(name, tree_dir))
        os.system("mv higgsCombine.{0}.Significance.mH120.root {1}/".format(name, tree_dir))
		        #run_combine = "combine -M Significance -n.{} {}/{}_nominal.txt -t -1 --expectSignal=1".format(name, input_dir, name)
		        #print run_combine
                #os.system(run_combine)
    
    for sample in samples:
        t = Thread(target=task, args=(sample,))
        t.start()
        #task(sample)
