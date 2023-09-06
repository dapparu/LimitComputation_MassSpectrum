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

    input_dir='datacards_SR3_test'
    tree_dir='limitTrees_SR3_test'
	
    os.system("mkdir -p {0}/".format(input_dir))
    os.system("mkdir -p {0}/".format(tree_dir))

    samples = [
        'Gluino800_2018',
        'Gluino1000_2018',
        'Gluino1400_2018',
        'Gluino1600_2018',
	    'Gluino1800_2018',
        'Gluino2000_2018',
        'Gluino2200_2018',
        'Gluino2400_2018',
        'Gluino2600_2018',
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
        run_combine = "combine -M AsymptoticLimits -n .{} -d {}/{}.txt --rRelAcc 0.000005 --rAbsAcc 0.000005 --rMin -1000.0 --rMax 1000.0".format(name, input_dir, name)
        print run_combine
        os.system(run_combine)
        run_combine = "combine -M Significance -n.{} {}/{}.txt".format(name, input_dir, name)
        print run_combine
        os.system(run_combine)
        os.system("mv higgsCombine.{0}.AsymptoticLimits.mH120.root {1}/".format(name, tree_dir))
        os.system("mv higgsCombine.{0}.Significance.mH120.root {1}/".format(name, tree_dir))
    
    for sample in samples:
        ### use these two lines to run on parallel
        t = Thread(target=task, args=(sample,))
        t.start()
        
        ### use this line to run one signal at once
        #task(sample)
