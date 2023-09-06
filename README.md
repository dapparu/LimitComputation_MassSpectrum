# LimitComputation_MassSpectrum
This repositery is dedicated to HSCP analysis, for the limits computation with the mass spectrum approach. 

## Setup working area

```bash
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_11_3_4
cd CMSSW_11_3_4/src/
cmsenv
```

For the following step you should have a ssh key associated to your GitHub account.
For more information, see [connecting-to-github-with-ssh-key](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent).

```bash
git clone -b master git@github.com:dapparu/git@github.com:dapparu/LimitComputation_MassSpectrum.git LimitComputation_MassSpectrum 
``` 

## Install the Combine packages and setup

```bash
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
```

Update to a recommended tag - currently the recommended tag is v9.1.0
More information on: https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/ 

```bash
cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v9.1.0
scramv1 b clean; scramv1 b # always make a clean build
```

Then install Combine Harvester package. 
More information on: https://cms-analysis.github.io/CombineHarvester/ 

```bash
cd CMSSW_11_3_4/src
cmsenv
git clone https://github.com/cms-analysis/CombineHarvester.git CombineHarvester
git checkout v2.0.0
scram b
```

Create ```HSCPLimit``` repository:

```bash
mkdir HSCPLimit
cd HSCPLimit
```
And finally install HSCP scripts: 

```bash 
git clone git@github.com:dapparu/LimitComputation_MassSpectrum.git LimitComputation_MassSpectrum
cd LimitComputation_MassSpectrum
```

## Import the cross-sections

The cross-section are already loaded in the directory ```xsec```. 
More information on: https://github.com/fuenfundachtzig/xsec

## Create the datacards

This part aims to create the datacards for each signal hypotheses. 

Use the script ```create_datacards.py```:

```bash
python create_datacards.py
```

Set the variable ```regionSignal``` to ```'SR1'```, ```'SR2'``` or ```'SR3'``` depending on what you want. 
Set the variables ```pathSignal``` and ```pathPred``` to open signal and predict mass distributions. The datacards are created for a whole bunch of different signal hypotheses and use the different predicted mass shapes due to systematics. 
Set the output directory with the variable ```outDataCardsDir```.
Setting the variable ```systSignal``` one will produce systematics budget for a given signal hypothesis, as a function of target masses (so systematics budget within the mass windows).

The bias correction parameters are currently hard-coded for each signal regions in both year. Be sure to update the values obtained with the package ```massSpectrum_bckgPrediction```, especially using the ```macroMass.py``` script (see: https://github.com/dapparu/massSpectrum_bckgPrediction). 

## Run Combine on the datacards

This part gives the way to run Combine limits and Combine significance on the previously produced datacards. 

Use the script ```run_datacards.py```:

```bash
python run_datacards.py
```

Set the ```input_dir``` and the ```tree_dir``` variables. The input directory of this code is the output directory of the previous stage.

By default, the code run the ```AsymptoticLimits``` and ```Significance``` methods, in parallel.

All the results are saved in the ```tree_dir``` directory.

## Produce the limits plots 

One uses the ```limit_plots.py``` script to produce limits plots:

```bash
python limit_plots.py
```

Set the ```labelSignal``` variable to extract limits on a given signal hypothesis. The supported hypothesis are the commented ones. 

The results are saved on the ```limit_plots_dir``` directory. 

## Significance computation

Work in progress... 

## Signal injection tests

Work in progress... 

## Mass shape analysis

Work in progress... 